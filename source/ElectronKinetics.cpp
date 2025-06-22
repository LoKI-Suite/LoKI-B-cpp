/** \file
 *
 *  Interfaces of classes that produce the EEDF and calculate swarm
 *  parameters and power terms.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   13 July 2023
 */

#include "LoKI-B/ElectronKinetics.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/EedfUtilities.h"
#include "LoKI-B/GridOps.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/OperatorsNew.h"
#include <cmath>
#include <limits>

//#define LOKIB_CREATE_SPARSITY_PICTURE
#ifdef LOKIB_CREATE_SPARSITY_PICTURE

#include "LoKI-B/Matrix2Picture.h"

#endif

// TODO [FUTURE]: Write a tridiagonal matrix class that stores the elements
//  in three separate vectors. It is desirable to do this in an Eigen compliant way, such
//  that these matrices can simply be added to dense matrices (the + and () operators are
//  the only operators to overload).

namespace loki
{

ElectronKinetics::ElectronKinetics(const std::filesystem::path &basePath, const json_type &cnf, WorkingConditions *workingConditions)
    : m_workingConditions(workingConditions),
    m_grid(Grid::fromConfig(cnf.at("numerics").at("energyGrid"))),
    mixture(basePath, &grid(), cnf, workingConditions),
    fieldOperator(grid()),
    inelasticOperator(grid()),
    carOperator(mixture.CARGases().empty() ? nullptr : new CAROperator(mixture.CARGases())),
    eedf(grid().nCells())
{
}

void ElectronKinetics::updateMaxEnergy(double uMax)
{
    m_grid.updateMaxEnergy(uMax);
}

void ElectronKinetics::updateMaxEnergyNonuniform(double uMax)
{
    m_grid.updateMaxEnergyNonuniform(uMax, mixture);
}

void ElectronKinetics::solve()
{
    doSolve();
}

ElectronKineticsBoltzmann::ElectronKineticsBoltzmann(const std::filesystem::path &basePath, const json_type &cnf, WorkingConditions *workingConditions)
: ElectronKinetics(basePath, cnf,workingConditions),
    ionizationOperator(getIonizationOperatorType(cnf.at("ionizationOperatorType"))),
    eeOperator(cnf.at("includeEECollisions") ? new ElectronElectronOperator(grid()) : nullptr),
    fieldMatrixSpatGrowth(grid().nCells(), grid().nCells()),
    ionSpatialGrowthD(grid().nCells(), grid().nCells()),
    ionSpatialGrowthU(grid().nCells(), grid().nCells()),
    alphaRedEff(0.),
    fieldMatrixTempGrowth(grid().nCells(), grid().nCells()),

    ionTemporalGrowth(grid().nCells(), grid().nCells()),
    CIEff(0.0),
    mixingParameter(cnf.at("numerics").at("nonLinearRoutines").at("mixingParameter")),
    growthModelType(getGrowthModelType(cnf.at("growthModelType"))),
    maxEedfRelError(cnf.at("numerics").at("nonLinearRoutines").at("maxEedfRelError")),
    maxPowerBalanceRelError(cnf.at("numerics").at("maxPowerBalanceRelError"))
{
    if (growthModelType == GrowthModelType::spatial && workingConditions->reducedExcFreqSI()!=0.0)
    {
        throw std::runtime_error("The excitation frequency must be zero when "
                                 "the spatial growth model is used.");
    }
    boltzmannMatrix.setZero(grid().nCells(), grid().nCells());

    /// \todo the following two tasks should probably be part of the IonizationOperator constructor
    ionizationOperator.ionConservativeMatrix.setZero(grid().nCells(), grid().nCells());
    if (ionizationOperator.ionizationOperatorType != IonizationOperatorType::conservative &&
        mixture.collision_data().hasCollisions(CollisionType::ionization))
    {
        ionizationOperator.ionizationMatrix.setZero(grid().nCells(), grid().nCells());
    }
    if (eeOperator)
    {
        eeOperator->initialize(grid());
    }

    // SPARSE INITIALIZATION
    /** \todo Could the matrices that are guaranteed to be diagonal just be Eigen::DiagonalMatrix?
     *        That saves a lot of time initializing and makes accessing easier (no 'coeffRef').
     *        Then the pattern-code below can be removed, only a resize needed. Question:
     *        does Eigen optimize matrix addition and multiplication for such matrices? That
     *        is not immediately clear to me.
     */
    std::vector<Eigen::Triplet<double>> diagPattern;
    diagPattern.reserve(grid().nCells());

    for (Grid::Index k = 0; k < grid().nCells(); ++k)
    {
        diagPattern.emplace_back(k, k, 0.);
    }
    std::vector<Eigen::Triplet<double>> tridiagPattern;
    tridiagPattern.reserve(3 * grid().nCells() - 2);

    for (Grid::Index k = 0; k < grid().nCells(); ++k)
    {
        if (k > 0)
            tridiagPattern.emplace_back(k - 1, k, 0.);

        tridiagPattern.emplace_back(k, k, 0.);

        if (k < grid().nCells() - 1)
            tridiagPattern.emplace_back(k + 1, k, 0.);
    }
    elasticMatrix.resize(grid().nCells(), grid().nCells());
    elasticMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());
    fieldMatrix.resize(grid().nCells(), grid().nCells()),
    fieldMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());
    if (carOperator.get())
    {
        CARMatrix.resize(grid().nCells(), grid().nCells());
        CARMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());
    }
    attachmentOperator.attachmentConservativeMatrix.resize(grid().nCells(), grid().nCells());
    attachmentOperator.attachmentConservativeMatrix.setZero(grid().nCells(), grid().nCells());
    attachmentOperator.attachmentMatrix.resize(grid().nCells(), grid().nCells());
    if (mixture.collision_data().hasCollisions(CollisionType::attachment))
    {
        attachmentOperator.attachmentMatrix.setFromTriplets(diagPattern.begin(), diagPattern.end());
    }
    if (growthModelType == GrowthModelType::spatial)
    {
        ionSpatialGrowthD.setFromTriplets(diagPattern.begin(), diagPattern.end());
        ionSpatialGrowthU.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());
        fieldMatrixSpatGrowth.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());
    }
    else if (growthModelType == GrowthModelType::temporal)
    {
        ionTemporalGrowth.setFromTriplets(diagPattern.begin(), diagPattern.end());
        fieldMatrixTempGrowth.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());
    }

    grid().updatedMaxEnergy.addListener(&ElectronKineticsBoltzmann::evaluateMatrix, this);
    m_workingConditions->updatedReducedField.addListener(&ElectronKineticsBoltzmann::evaluateFieldOperator, this);
    /// \todo Do this here? Not all parameters may have been set at this point.
    this->evaluateMatrix();
}

void ElectronKineticsBoltzmann::evaluateFieldOperator()
{
    const double WoN = m_workingConditions->reducedExcFreqSI();
    /// \todo Should we use the real value (CIEff) here? That will be 0 for DC or spatial growth
    const double dummyCIEff = 0.0;
    fieldOperator.evaluate(grid(),mixture.collision_data().totalCrossSection(),WoN,dummyCIEff,fieldMatrix);
}

void ElectronKineticsBoltzmann::evaluateMatrix()
{
    const double Tg = m_workingConditions->gasTemperature();
    mixture.collision_data().evaluateTotalAndElasticCS(grid());

    elasticOperator.evaluate(grid(),mixture.collision_data().elasticCrossSection(),Tg,elasticMatrix);

    evaluateFieldOperator();

    // if (carOperator.get())
    // {
    //     carOperator->evaluate(grid(),Tg,CARMatrix);
    // }

    inelasticOperator.evaluateInelasticOperators(grid(),mixture);

    // if (mixture.collision_data().hasCollisions(CollisionType::ionization))
    //     ionizationOperator.evaluateIonizationOperator(grid(),mixture);

    // if (mixture.collision_data().hasCollisions(CollisionType::attachment))
    //     attachmentOperator.evaluateAttachmentOperator(grid(),mixture);

    // if (eeOperator)
    // {
    //     eeOperator->initialize(grid());
    // }
}

void ElectronKineticsBoltzmann::doSolve()
{
    if (grid().smartGrid())
    {
//#define LOKIB_USE_BISECTING_SMART_GRID
#ifdef LOKIB_USE_BISECTING_SMART_GRID
        solveSmartGrid2();
#else
        solveSmartGrid();
#endif
    }
    else
    {
        solveSingle();
    }
    evaluatePower();
    // we finished the solution procedure. If the power balance is
    // still not good, issue a warning.
    if (power.relativeBalance > maxPowerBalanceRelError)
    {
        Log<PowerBalanceError>::Warning(maxPowerBalanceRelError);
    }

    // evaluate derived data
    mixture.collision_data().evaluateRateCoefficients(grid(),eedf);
    evaluateSwarmParameters();
    evaluateFirstAnisotropy();
    obtainedNewEedf.emit(grid(), eedf, *m_workingConditions, power, mixture.collision_data(), swarmParameters,
                         &firstAnisotropy);
#ifdef LOKIB_CREATE_SPARSITY_PICTURE
    const std::string xpm_fname{"system_matrix.xpm"};
    std::cout << "Creating '" << xpm_fname << "'." << std::endl;
    writeXPM(boltzmannMatrix, xpm_fname);
#endif
}

void ElectronKineticsBoltzmann::invertLinearMatrixNew()
{
    experimental::ElasticOperator elastic_operator(grid());
    experimental::FieldOperator field_operator(grid());

    elastic_operator.evaluate(
        grid(),
        mixture.collision_data().elasticCrossSection(),
        m_workingConditions->gasTemperature()
    );
    field_operator.evaluate(
        grid(),
        mixture.collision_data().totalCrossSection(),
        m_workingConditions->reducedFieldSI()
    );

    const auto drift_coeff = elastic_operator.drift_coefficient() + field_operator.drift_coefficient();
    const auto diff_coeff = elastic_operator.diffusion_coefficient() + field_operator.diffusion_coefficient();

    Vector peclet = drift_coeff.array() * grid().duNodes().array() / diff_coeff.array();

    Matrix baseMatrix(grid().nCells(), grid().nCells());

    baseMatrix.setZero();

    // Apply the Scharfetter-Gummel discretization scheme with a zero-flux
    // boundary condition on both boundaries.
    for (Grid::Index i = 0; i < grid().nCells(); i++) {
        if (i > 0) {
            baseMatrix(i, i) -= drift_coeff[i] / (1. - std::exp(peclet[i]));
            baseMatrix(i, i - 1) -= drift_coeff[i] / (1. - std::exp(-peclet[i]));
        }

        if (i < grid().nCells() - 1) {
            baseMatrix(i, i) += drift_coeff[i + 1] / (1. - std::exp(-peclet[i + 1]));
            baseMatrix(i, i + 1) += drift_coeff[i + 1] / (1. - std::exp(peclet[i + 1]));
        }
    }
    // Central difference scheme.
    // for (Grid::Index i = 0; i < grid().nCells(); i++) {
    //     if (i > 0) {
    //         baseMatrix(i, i) -= drift_coeff[i] / 2. - diff_coeff[i] / grid().duNode(i);
    //         baseMatrix(i, i - 1) -= drift_coeff[i] / 2. + diff_coeff[i] / grid().duNode(i);
    //     }

    //     if (i < grid().nCells() - 1) {
    //         baseMatrix(i, i) += drift_coeff[i + 1] / 2. + diff_coeff[i + 1] / grid().duNode(i + 1);
    //         baseMatrix(i, i + 1) += drift_coeff[i + 1] / 2. - diff_coeff[i + 1] / grid().duNode(i + 1);
    //     }
    // }

    // TODO: Determine a suitable value for the electron temperature.
    eedf = makePrescribedEDF(grid(), 1, 0.1, false);

    Vector eedf_cur(eedf);

    double error = std::numeric_limits<double>::max();

    for (int i = 0; i < 100 && error > 1e-5; ++i) {
        eedf_cur = eedf;

        inelasticOperator.evaluateInelasticOperatorsNew(grid(), eedf, this->mixture);
        boltzmannMatrix = baseMatrix + inelasticOperator.inelasticMatrix;

        // Apply the constant logarithmic decay boundary condition.
        // const auto N = grid().nCells() - 1;
        // boltzmannMatrix.row(N).setZero();
        // boltzmannMatrix(N, N) = 1;
        // boltzmannMatrix(N, N - 1) = -std::pow(eedf[N - 1] / eedf[N - 2], grid().duNode(N - 1) / grid().duNode(N - 2));

        invertMatrix(boltzmannMatrix);

        error = (eedf - eedf_cur).cwiseQuotient(eedf).cwiseAbs().sum();
        Log<Message>::Warning("EEDF rel error: ", error);
    }

    normalizeEDF(eedf, grid());
}

void ElectronKineticsBoltzmann::invertLinearMatrix()
{
    const double EoN = m_workingConditions->reducedFieldSI();
    // Here the Conservative ionization and attachment matrices are added.
    // Otherwise, the sum is the same as in solveSpatialGrowthMatrix
    // (note that in solveTemporalGrowthMatrix fieldMatrix is not added).
    boltzmannMatrix
        = elasticMatrix
        + fieldMatrix*(EoN*EoN)
        + inelasticOperator.inelasticMatrix;
    //     + ionizationOperator.ionConservativeMatrix
    //     + attachmentOperator.attachmentConservativeMatrix;
    // if (carOperator)
    // {
    //     boltzmannMatrix += CARMatrix;
    // }

    invertMatrix(boltzmannMatrix);
    normalizeEDF(eedf, grid());
}

void ElectronKineticsBoltzmann::invertMatrix(Matrix &matrix)
{
    solveEEDF(eedf,matrix,grid());
}

void ElectronKineticsBoltzmann::solveSpatialGrowthMatrix()
{
    const double EoN = m_workingConditions->reducedFieldSI();

    const Vector& cellTotalCrossSection(mixture.collision_data().totalCellCrossSection());
    boltzmannMatrix
        = elasticMatrix
        + fieldMatrix*(EoN*EoN)
        + inelasticOperator.inelasticMatrix
        + ionizationOperator.ionizationMatrix
        + attachmentOperator.attachmentMatrix;
    if (carOperator)
    {
        boltzmannMatrix += CARMatrix;
    }
    // *add* the discretization of the ee term
    if (eeOperator)
    {
        eeOperator->discretizeTerm(boltzmannMatrix,grid());
    }

    /* Store the diagonals of the Boltzmann matrix that are going to be
     * affected by the growth terms, so we can restore these later,
     * before updated growth terms are going to be applied.
     */
    Vector baseDiag(grid().nCells()), baseSubDiag(grid().nCells()), baseSupDiag(grid().nCells());

    for (Grid::Index k = 0; k < grid().nCells(); ++k)
    {
        baseDiag[k] = boltzmannMatrix(k, k);

        if (k > 0)
            baseSubDiag[k] = boltzmannMatrix(k, k - 1);

        if (k < grid().nCells() - 1)
            baseSupDiag[k] = boltzmannMatrix(k, k + 1);
    }

    const Vector coefsCI = SI::gamma * grid().duCells().transpose()
                    * (ionizationOperator.ionizationMatrix + attachmentOperator.attachmentMatrix);

    double CIEffNew = eedf.dot(coefsCI);

    /** \todo DB: diffusion and mobility components of the spatial growth terms
     *  can be removed since this is already done in the directMixing function
     *  Jvd: Daan, can you explain this statement? What can be done/simplified?
     *       These matrices are set up, assigned values to and used in the
     *       spatial growth case.
     */
    ionSpatialGrowthD.setZero();
    ionSpatialGrowthU.setZero();

    /* In the following block we evaluate D_e and mu_eE as given by equations
     * 19a and 19b of \cite Tejero2019. Note that we are handling the SST case
     * in this function, and for x=SST we have Omega_x=sigma_c (equation 5b).
     * We first define D0(u) = (1./3)*u/Omega_SST(u) = (1./3)*u/sigma_c(u).
     *
     * Subsequently, we evaluate D_eN=+gamma*int_0^infty D0(u)f(u)du (eqn. 19a)
     * and mu_eE=-gamma*int_0^infty D0(u)[df/du]du (19b).
     */
    const Vector D0 = grid().getCells().array() / (3. * cellTotalCrossSection).array();
    const Vector D0Nodes = grid().getNodes().array() / (3. * mixture.collision_data().totalCrossSection()).array();

    // This is 33a from \cite Manual_1_0_0
    double ND  =   SI::gamma* energyIntegral(grid(),D0,eedf);
    double muE = - SI::gamma * fNodegPrimeEnergyIntegral(grid(), D0Nodes, eedf) * EoN;

    double alphaRedEffOld = 0.;

    /* The initial guess for the eedf may lead a negative discriminant.
     * In this case the reduced townsend coefficient should be calculated
     * on the basis of the assumption that there is no electron density
     * gradient.
     */
    /** \todo discr==0 corresponds to muE/2ND = 2*CIEffNew/muE
     *        and this will be assigned to alphaRedEffNew by the discr>=0
     *        code path. But for discr<0 we assign CIEffNew / muE, resulting
     *        in a discontinuity. It is almost as if a factor 2 is missing
     *        somewhere. (NB: this can in theory frustrate convergence
     *        for the case that discr is around 0, I think.
     */
    /** \todo The following line implements the solution of equation 22 of
     *  \cite Tejero2019, but there are two differences:
     *   - the name alpha_eff is used in the paper, whereas the word *reduced*
     *     townsend coefficient is used in the code. Better be consistent.
     *   - Next, the case that the roots are complex is not discussed in the
     *     text, only here in the code. It would be nice to understand a bit
     *     better underwhat circumstances that happens, and more importantly:
     *     why that does *not* happen for a 'correct' f(u).
     */
    const double discriminant = muE*muE - 4*CIEffNew*ND;
    double alphaRedEffNew = (discriminant < 0.)
        ? CIEffNew / muE
        : (muE - std::sqrt(discriminant)) / (2 * ND);

    uint32_t iter = 0;
    bool hasConverged = false;

    // note: this g_fieldSpatialBase = (E/N)*D^0(u)
    const Vector g_fieldSpatialBase = (EoN / 3) * grid().getNodes().array() / mixture.collision_data().totalCrossSection().array();

    while (!hasConverged)
    {
        // note: this g_fieldSpatialGrowth = (alphaEffNew/N)*(E/N)*D^0(u)
        g_fieldSpatialGrowth = alphaRedEffNew * g_fieldSpatialBase;
        g_fieldSpatialGrowth[0] = 0.;
        g_fieldSpatialGrowth[grid().nCells()] = 0.;

        /* note: what is still missing is (alpha/N)*(E/N)*D0*df/du := C_k*df/du,
         * which is represented by ionSpatialGrowthU.
         * In an *internal* point k we have C_(df/du)_k \approx C_k(f_{k+1}-f_{k-1}})/(2*du),
         * with C_k = (alpha/N)*(E/N)*D0_k. This expression is used when
         * USE_D0_FOR_ionSpatialGrowthU is defined to 1. See the todo below
         * for a note on the boundary cells.
         *
         * NOTE:
         * In the MATLAB code, +/-C_k is expressed in terms of U0{sup,inf},
         * which makes the expression a bit more difficult to understand.
         * What is the reason for that? Consistency with the evaluation of
         * mu_eE? That can also be achieved without using U{inf,Usup}.
         * To understand the original expressions below, note that:
         *   alphaRedEffNew * U0inf[k-1] == alphaRedEffNew*(-EoN / (2. * grid().du()) * D0[k]) = -C_k/(2*du)
         *   alphaRedEffNew * U0sup[k+1] == alphaRedEffNew*(+EoN / (2. * grid().du()) * D0[k]) = +C_k/(2*du),
         */
        /** \todo At boundary points, we do not seem to be implementing the term correctly.
         * The problem (or misunderstanding on my side (JvD)) is similar to
         * that in the evaluation of mu_eE, see elsewhere.
         */
        if (grid().isUniform())
        {
            for (Grid::Index k = 0; k < grid().nCells(); ++k)
            {
                boltzmannMatrix(k, k) = baseDiag[k];

                if (k > 0)
                    boltzmannMatrix(k,k-1) = baseSubDiag[k];

                if (k < grid().nCells() - 1)
                    boltzmannMatrix(k,k+1) = baseSupDiag[k];
            }

            for (Grid::Index k = 0; k < grid().nCells(); ++k)
            {
                /* Handle ionSpatialGrowthD, which is defined such that
                * [ionSpatialGrowthD*eedf]_k = (alphaEffNew/N)^2*[D0*f]_k.
                * (This only has a diagonal element.)
                */
                ionSpatialGrowthD.coeffRef(k, k) = alphaRedEffNew * alphaRedEffNew * D0[k];
                boltzmannMatrix(k, k) += ionSpatialGrowthD.coeff(k, k);

                /* Handle ionSpatialGrowthU, which is defined such that
                * [ionSpatialGrowthU*eedf]_k = (alpha/N)*(E/N)*[D0*df/du]_k.
                */
                if (k==0)
                {
                    ionSpatialGrowthU.coeffRef(k, k    ) = -alphaRedEffNew*EoN*D0Nodes[k+1] / (2.*grid().du());
                    ionSpatialGrowthU.coeffRef(k, k + 1) = +alphaRedEffNew*EoN*D0Nodes[k+1] / (2.*grid().du());
                    boltzmannMatrix(k,k  ) += ionSpatialGrowthU.coeff(k,k);
                    boltzmannMatrix(k,k+1) += ionSpatialGrowthU.coeff(k,k+1);
                }
                else if (k==grid().nCells() - 1)
                {
                    ionSpatialGrowthU.coeffRef(k, k - 1) = -alphaRedEffNew*EoN*(D0Nodes[k] + D0Nodes[k+1]) / (2. * grid().du());
                    ionSpatialGrowthU.coeffRef(k, k    ) = +alphaRedEffNew*EoN*(D0Nodes[k] + D0Nodes[k+1]) / (2. * grid().du());
                    boltzmannMatrix(k,k-1) += ionSpatialGrowthU.coeff(k,k-1);
                    boltzmannMatrix(k,k  ) += ionSpatialGrowthU.coeff(k,k);
                }
                else
                {
                    ionSpatialGrowthU.coeffRef(k, k) = -alphaRedEffNew*EoN*(D0Nodes[k+1] - D0Nodes[k]) / (2.*grid().du());
                    ionSpatialGrowthU.coeffRef(k, k - 1) = -alphaRedEffNew*EoN*D0Nodes[k] / (2.*grid().du());
                    ionSpatialGrowthU.coeffRef(k, k + 1) = +alphaRedEffNew*EoN*D0Nodes[k+1] / (2.*grid().du());
                    boltzmannMatrix(k,k  ) += ionSpatialGrowthU.coeff(k,k  );
                    boltzmannMatrix(k,k-1) += ionSpatialGrowthU.coeff(k,k-1);
                    boltzmannMatrix(k,k+1) += ionSpatialGrowthU.coeff(k,k+1);
                }

                /* Handle fieldMatrixSpatGrowth, which is defined such that
                * [fieldMatrixSpatGrowth*eedf]_k = (alphaEffNew/N)*(E/N)*[d(D^0*f0)/du]_k.
                */
                fieldMatrixSpatGrowth.coeffRef(k, k) = (g_fieldSpatialGrowth[k + 1] - g_fieldSpatialGrowth[k]) / (2*grid().du());
                boltzmannMatrix(k, k) += fieldMatrixSpatGrowth.coeff(k, k);
                if (k > 0)
                {
                    fieldMatrixSpatGrowth.coeffRef(k, k - 1) = -g_fieldSpatialGrowth[k] / (2*grid().du());
                    boltzmannMatrix(k, k - 1) += fieldMatrixSpatGrowth.coeff(k, k - 1);
                }

                if (k < grid().nCells() - 1)
                {
                    fieldMatrixSpatGrowth.coeffRef(k, k + 1) = g_fieldSpatialGrowth[k + 1] / (2*grid().du());
                    boltzmannMatrix(k, k + 1) += fieldMatrixSpatGrowth.coeff(k, k + 1);
                }
            }
        } else
        {
            for (Grid::Index k = 0; k < grid().nCells(); ++k)
            {
                boltzmannMatrix(k, k) = baseDiag[k];

                if (k > 0)
                    boltzmannMatrix(k,k-1) = baseSubDiag[k];

                if (k < grid().nCells() - 1)
                    boltzmannMatrix(k,k+1) = baseSupDiag[k];
            }

            for (Grid::Index k = 0; k < grid().nCells(); ++k)
            {
                /* Handle ionSpatialGrowthD, which is defined such that
                * [ionSpatialGrowthD*eedf]_k = (alphaEffNew/N)^2*[D0*f]_k.
                * (This only has a diagonal element.)
                */
                ionSpatialGrowthD.coeffRef(k, k) = alphaRedEffNew * alphaRedEffNew * D0[k];
                boltzmannMatrix(k, k) += ionSpatialGrowthD.coeff(k, k);

                /* Handle ionSpatialGrowthU, which is defined such that
                * [ionSpatialGrowthU*eedf]_k = (alpha/N)*(E/N)*[D0*df/du]_k.
                */
                if (k==0)
                {
                    ionSpatialGrowthU.coeffRef(k, k    ) = -alphaRedEffNew*EoN*D0[k] / (grid().duNode(1));
                    ionSpatialGrowthU.coeffRef(k, k + 1) = +alphaRedEffNew*EoN*D0[k] / (grid().duNode(1));
                    boltzmannMatrix(k,k  ) += ionSpatialGrowthU.coeff(k,k);
                    boltzmannMatrix(k,k+1) += ionSpatialGrowthU.coeff(k,k+1);
                }
                else if (k==grid().nCells() - 1)
                {
                    ionSpatialGrowthU.coeffRef(k, k - 1) = -alphaRedEffNew*EoN*D0[k] / (grid().duNode(k));
                    ionSpatialGrowthU.coeffRef(k, k    ) = +alphaRedEffNew*EoN*D0[k] / (grid().duNode(k));
                    boltzmannMatrix(k,k-1) += ionSpatialGrowthU.coeff(k,k-1);
                    boltzmannMatrix(k,k  ) += ionSpatialGrowthU.coeff(k,k);
                }
                else
                {
                    ionSpatialGrowthU.coeffRef(k, k - 1) = -alphaRedEffNew*EoN*D0[k] / (grid().duNode(k) + grid().duNode(k + 1));
                    ionSpatialGrowthU.coeffRef(k, k + 1) = +alphaRedEffNew*EoN*D0[k] / (grid().duNode(k) + grid().duNode(k + 1));
                    boltzmannMatrix(k,k-1) += ionSpatialGrowthU.coeff(k,k-1);
                    boltzmannMatrix(k,k+1) += ionSpatialGrowthU.coeff(k,k+1);
                }

                /* Handle fieldMatrixSpatGrowth, which is defined such that
                * [fieldMatrixSpatGrowth*eedf]_k = (alphaEffNew/N)*(E/N)*[d(D^0*f0)/du]_k.
                */
                fieldMatrixSpatGrowth.coeffRef(k, k) = (g_fieldSpatialGrowth[k+1] - g_fieldSpatialGrowth[k]) / grid().duCell(k);
                boltzmannMatrix(k, k) += fieldMatrixSpatGrowth.coeff(k, k);

                if (k > 0) {
                    fieldMatrixSpatGrowth.coeffRef(k, k-1) =
                        - (g_fieldSpatialGrowth[k] + g_fieldSpatialGrowth[k+1])
                        / 2.
                        / (grid().duNode(k) + grid().duNode(k+1));
                    boltzmannMatrix(k, k-1) += fieldMatrixSpatGrowth.coeff(k, k-1);
                }

                if (k < grid().nCells() - 1) {
                    fieldMatrixSpatGrowth.coeffRef(k, k+1) =
                        (g_fieldSpatialGrowth[k] + g_fieldSpatialGrowth[k+1])
                        / 2.
                        / (grid().duNode(k) + grid().duNode(k+1));
                    boltzmannMatrix(k, k+1) += fieldMatrixSpatGrowth.coeff(k, k+1);
                }
            }
        }
        Vector eedfNew = eedf;

        invertMatrix(boltzmannMatrix);

        CIEffNew = eedf.dot(coefsCI);

        ND = SI::gamma * energyIntegral(grid(), D0, eedf);
        muE = -SI::gamma * fNodegPrimeEnergyIntegral(grid(), D0Nodes, eedf) * EoN;

#if 0
        std::cout << "muE: " << muE << std::endl;
        std::cout << "alt: " << SI::gamma*grid().duCells().dot(grid().getCells().asDiagonal()*(fieldMatrix*eedf))*EoN << std::endl;
#endif

        /** \todo See the notes just above this loop for a note about the discontinuity.
         */
        const double discriminant_new = muE * muE - 4 * CIEffNew * ND;

        alphaRedEffOld = alphaRedEffNew;
        alphaRedEffNew = (discriminant_new < 0) ? CIEffNew / muE : (muE - std::sqrt(discriminant_new)) / (2 * ND);

        if (((alphaRedEffNew == 0 || std::abs(alphaRedEffNew - alphaRedEffOld) / alphaRedEffOld < 1.e-9) &&
             maxRelDiff(eedfNew,eedf) < maxEedfRelError) ||
            iter > 150)
        {
            hasConverged = true;

            /** There is no maximum number of iterations if ee collisions are enabled.
             *  Is there a reason for that? This is not very safe, it seems.
             */
            if (iter > 150 && !eeOperator)
                Log<Message>::Warning("Iterative spatial growth scheme did not converge.");
        }

        alphaRedEffNew = mixingParameter * alphaRedEffNew + (1 - mixingParameter) * alphaRedEffOld;

        ++iter;
    }

    std::cerr << "Spatial growth routine converged in: " << iter << " iterations.\n";

    alphaRedEff = alphaRedEffOld;
}

void ElectronKineticsBoltzmann::solveTemporalGrowthMatrix()
{
    // same sum as in solveSpatialGrowthMatrix, except that fieldMatrix is skipped
    boltzmannMatrix
        = elasticMatrix
        + inelasticOperator.inelasticMatrix
        + ionizationOperator.ionizationMatrix
        + attachmentOperator.attachmentMatrix;
    if (carOperator)
    {
        boltzmannMatrix += CARMatrix;
    }
    // *add* the discretization of the ee term
    if (eeOperator)
    {
        eeOperator->discretizeTerm(boltzmannMatrix,grid());
    }

    const double EoN = m_workingConditions->reducedFieldSI();
    const double WoN = m_workingConditions->reducedExcFreqSI();

    Vector baseDiag(grid().nCells()), baseSubDiag(grid().nCells()), baseSupDiag(grid().nCells());

    for (Grid::Index k = 0; k < grid().nCells(); ++k)
    {
        baseDiag[k] = boltzmannMatrix(k, k);

        if (k > 0)
            baseSubDiag[k] = boltzmannMatrix(k, k - 1);

        if (k < grid().nCells() - 1)
            baseSupDiag[k] = boltzmannMatrix(k, k + 1);
    }

    // CIEff is <nu_eff>/N
    const Vector coefsCI = SI::gamma * grid().duCells().transpose()
                    * (ionizationOperator.ionizationMatrix + attachmentOperator.attachmentMatrix);
    double CIEffNew = eedf.dot(coefsCI);
    double CIEffOld = 0;

    Vector eedfNew(grid().nCells());

    bool hasConverged = false;
    uint32_t iter = 0;

    while (!hasConverged)
    {
 //       Log<Message>::Notify("Iteration ", iter);

        // CIEff is <nu_eff>/N, so growthFactor = <nu_eff>/(N*gamma)
        fieldOperator.evaluate(grid(),mixture.collision_data().totalCrossSection(),WoN,CIEffNew,fieldMatrixTempGrowth);
        const long double growthFactor = CIEffNew / SI::gamma;
        for (Grid::Index k = 0; k < grid().nCells(); ++k)
        {
            // Manual 2.2.0, 7a (with an minus sign because all terms are negated):
            ionTemporalGrowth.coeffRef(k, k) = -growthFactor * std::sqrt(grid().getCell(k));
            boltzmannMatrix(k, k) = baseDiag[k] + fieldMatrixTempGrowth.coeff(k, k)*EoN*EoN +
                                                           ionTemporalGrowth.coeff(k, k);

            if (k > 0)
            {
                boltzmannMatrix(k, k - 1) = baseSubDiag[k] + fieldMatrixTempGrowth.coeff(k, k - 1)*EoN*EoN;
            }

            if (k < grid().nCells() - 1)
            {
                boltzmannMatrix(k, k + 1) = baseSupDiag[k] + fieldMatrixTempGrowth.coeff(k, k + 1)*EoN*EoN;
            }

        }

        eedfNew = eedf;

        invertMatrix(boltzmannMatrix);

        CIEffOld = CIEffNew;
        CIEffNew = eedf.dot(coefsCI);

        if (((CIEffNew == 0 || std::abs(CIEffNew - CIEffOld) / CIEffOld < 1e-11) &&
             maxRelDiff(eedfNew,eedf) < maxEedfRelError) ||
            iter > 150)
        {
            hasConverged = true;

            if (iter > 150 && !eeOperator)
                Log<Message>::Warning("Iterative temporal growth scheme did not converge.");
        }

        CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

        ++iter;
    }

 //   std::cerr << "Temporal growth routine converged in: " << iter << " iterations.\n";

    CIEff = CIEffOld;
}

void ElectronKineticsBoltzmann::solveEEColl()
{
    assert(eeOperator);
    // Splitting all possible options for the best performance.

    const double EoN = m_workingConditions->reducedFieldSI();
    /** \todo What if only one of the following is true? Then we use e.g. the ionizationMatrix, while includeNonConservativeIonization==false.
     *  Is that correct?
     */
    if (ionizationOperator.includeNonConservativeIonization || attachmentOperator.includeNonConservativeAttachment)
    {
        /** \todo For spatial growth, we add fieldMatrix AND fieldMatrixSpatGrowth,
         *  for temporal growth only fieldMatrixTempGrowth. Is that correct?
         */
        if (growthModelType == GrowthModelType::spatial)
        {
            boltzmannMatrix
                = ionizationOperator.ionizationMatrix
                + attachmentOperator.attachmentMatrix
                + elasticMatrix
                + inelasticOperator.inelasticMatrix
                + fieldMatrix*(EoN*EoN)
                + ionSpatialGrowthD
                + ionSpatialGrowthU
                + fieldMatrixSpatGrowth;
        }
        else if (growthModelType == GrowthModelType::temporal)
        {
            boltzmannMatrix
                = ionizationOperator.ionizationMatrix
                + attachmentOperator.attachmentMatrix
                + elasticMatrix
                + inelasticOperator.inelasticMatrix
                + ionTemporalGrowth
                + fieldMatrixTempGrowth*(EoN*EoN);
        }
    }
    else
    {
        boltzmannMatrix
                = ionizationOperator.ionConservativeMatrix
                + attachmentOperator.attachmentConservativeMatrix
                + elasticMatrix
                + inelasticOperator.inelasticMatrix
                + fieldMatrix*(EoN*EoN);
    }
    if (carOperator)
    {
        boltzmannMatrix += CARMatrix;
    }

    /// \todo Make the remainder of this function a member of the eeOperator, passing the boltzmannMatrix as arument?

    const double ne = m_workingConditions->electronDensity();
    const double n0 = m_workingConditions->gasDensity();

    // Storing the initial diagonals of the matrix in three separate vectors.
    // This allows us to skip the usage of 'matrixAux', saving a good amount of
    // memory for bigger matrices.

    Vector baseDiag(grid().nCells()), baseSubDiag(grid().nCells()), baseSupDiag(grid().nCells());

    for (Grid::Index k = 0; k < grid().nCells(); ++k)
    {
        baseDiag[k] = boltzmannMatrix(k, k);

        if (k > 0)
            baseSubDiag[k] = boltzmannMatrix(k, k - 1);

        if (k < grid().nCells() - 1)
            baseSupDiag[k] = boltzmannMatrix(k, k + 1);
    }

    double ratioNew = 0.;
    Vector eedfNew = eedf;

    bool hasConverged = false;
    uint32_t iter = 0;

    while (!hasConverged)
    {
        eeOperator->update_g_ee_AB(grid(),eedf,ne,n0);
        // restore boltzmannMatrix to the situation without ee collisions
        for (Grid::Index k = 0; k < grid().nCells(); ++k)
        {
            boltzmannMatrix(k, k) = baseDiag[k];

            if (k > 0)
                boltzmannMatrix(k, k - 1) = baseSubDiag[k];

            if (k < grid().nCells() - 1)
                boltzmannMatrix(k, k + 1) = baseSupDiag[k];
        }
        // *add* the discretization of the ee term
        eeOperator->discretizeTerm(boltzmannMatrix,grid());

        invertMatrix(boltzmannMatrix);

        evaluatePower();

        const double ratio = std::abs(power.electronElectronNet / power.reference);

        if (maxRelDiff(eedfNew,eedf) < maxEedfRelError)
        {
            if (std::abs(ratio) < 1.e-9)
            {
                hasConverged = true;
            }
            else if (std::abs(ratio) > 1.e-9 && iter > 200)
            {
                Log<PowerRatioError>::Warning(std::abs(ratio));
                hasConverged = true;
            }
        }
        else if (iter == 300 && !(attachmentOperator.includeNonConservativeAttachment || ionizationOperator.includeNonConservativeIonization))
        {
            hasConverged = true;
            Log<Message>::Warning("Electron-electron iterative scheme did not converge.");
        }

        if (iter > 0 && !hasConverged)
        {
            double ratioOld = ratioNew;
            ratioNew = ratio;

            const Vector eedfOld = eedfNew;
            eedfNew = eedf;

            eedf = eedfNew - (ratioNew / (ratioNew - ratioOld)) * (eedfNew - eedfOld);

            /** \todo Is this an appropriate clipping criterium?
             *        Why take the absolute value, not some fixed small number?
             */
            for (Grid::Index i = 0; i < grid().nCells(); ++i)
            {
                if (eedf[i] < 0)
                {
                    eedf[i] = std::abs(eedf[i]);
                }
            }
        }
        iter++;
    }

    std::cerr << "e-e routine converged in: " << iter << " iterations.\n";
}

void ElectronKineticsBoltzmann::obtainTimeIndependentSolution()
{
    // invertLinearMatrix();
    invertLinearMatrixNew();
    return;

    /* Maybe we are done. But we need to do more work if non-linear terms are
     * present:
     *   1) growth terms (with or without ee collisions)
     *   2) ee collisions.
     * these cases are handled below.
     */
    const bool includeGrowthModel = attachmentOperator.includeNonConservativeAttachment || ionizationOperator.includeNonConservativeIonization;
    if (includeGrowthModel)
    {
        void (ElectronKineticsBoltzmann::*growthFunc)() = nullptr;
        switch (growthModelType)
        {
        case GrowthModelType::spatial:
            ionSpatialGrowthD.setZero();
            ionSpatialGrowthU.setZero();
            fieldMatrixSpatGrowth.setZero();

            growthFunc = &ElectronKineticsBoltzmann::solveSpatialGrowthMatrix;
            break;

        case GrowthModelType::temporal:
            ionTemporalGrowth.setZero();
            fieldMatrixTempGrowth.setZero();

            growthFunc = &ElectronKineticsBoltzmann::solveTemporalGrowthMatrix;
            break;
        }
        if (eeOperator)
        {
            /* This is the situation that is visualized in the flow chart (figure 3)
             * in the manual. The calculation without non-linear terms has already
             * been done at the start of this function. We now enter a double loop.
             * Inside the main loop 1) the growth model is run and 2) the ee model
             * is run. The main loop finishes when the last ee-run no longer changes
             * the EEDF (beyond tolerance settings).
             */

            /* We start without ee collisions. The eeOperator is updated inside the
             * call to solveEEColl(); the updated eeOperator will be used in the
             * next call to the growth function.
             */
            eeOperator->clear();
            Vector eedfOld;
            uint32_t globalIter = 0;
            /** \todo Should this be configurable? Note that this is also
             *  hardcoded to 20 in MATLAB
             */
            const uint32_t maxGlobalIter = 20;
            while (globalIter < maxGlobalIter)
            {
                (this->*growthFunc)();
                eedfOld = eedf;

                solveEEColl();

                if (maxRelDiff(eedf,eedfOld) < maxEedfRelError)
                    break;

                ++globalIter;

                if (globalIter == maxGlobalIter)
                    Log<GlobalIterError>::Warning(globalIter);
            }
        }
        else
        {
            (this->*growthFunc)();
        }
    }
    else if (eeOperator)
    {
            Log<Message>::Notify("Starting e-e collision routine.");
            solveEEColl();
    }
#if 0
    evaluatePower();
    const double EoN = m_workingConditions->reducedFieldSI();
    const double elP = SI::gamma*grid().duCells().dot(grid().getCells().asDiagonal()*(elasticMatrix*eedf));
    const double efP = SI::gamma*grid().duCells().dot(grid().getCells().asDiagonal()*(fieldMatrix*eedf))*(EoN*EoN);
    const double efPSG = SI::gamma*grid().duCells().dot(grid().getCells().asDiagonal()*(fieldMatrixSpatGrowth*eedf));
    std::cout << "Prat: " << elP/power.elasticNet << ", Prat_field = " << (efP+efPSG)/power.field << std::endl;
    std::cout << "P: " << efP << ", sg = " << efPSG << ", sum = " << power.field << std::endl;
#endif
}

void ElectronKineticsBoltzmann::solveSingle()
{
    obtainTimeIndependentSolution();
}

void ElectronKineticsBoltzmann::solveSmartGrid()
{
    assert (grid().smartGrid());
    const Grid::SmartGridParameters& smartGrid = *grid().smartGrid();
    solveSingle();
    double decades = calcDecades(eedf[0],eedf[grid().nCells()-1]);
    //std::cout << "decades: " << decades << ", uMax: " << grid().uMax() << std::endl;

    while (decades < smartGrid.minEedfDecay)
    {
        if (grid().isUniform())
        {
            updateMaxEnergy(grid().uMax() * (1 + smartGrid.updateFactor));
        } else
        {
            updateMaxEnergyNonuniform(grid().uMax() * (1 + smartGrid.updateFactor));
        }
        solveSingle();
        decades = calcDecades(eedf[0],eedf[grid().nCells()-1]);
        //std::cout << "decades: " << decades << ", uMax: " << grid().uMax() << std::endl;
    }

    while (decades > smartGrid.maxEedfDecay)
    {
        if (grid().isUniform())
        {
            updateMaxEnergy(grid().uMax() / (1 + smartGrid.updateFactor));
        } else
        {
            updateMaxEnergyNonuniform(grid().uMax() / (1 + smartGrid.updateFactor));
        }
        solveSingle();
        decades = calcDecades(eedf[0],eedf[grid().nCells()-1]);
        //std::cout << "decades: " << decades << ", uMax: " << grid().uMax() << std::endl;
    }
}

void ElectronKineticsBoltzmann::solveSmartGrid2()
{
    assert (grid().smartGrid());
    const Grid::SmartGridParameters& smartGrid = *grid().smartGrid();

    // 0. Set uM and uP equal to the present (initial) uMax, solve
    //    and calculate decades.
    double uM=grid().uMax();
    double uP=grid().uMax();
    solveSingle();
    double decades = calcDecades(eedf[0],eedf[grid().nCells()-1]);
    std::cout << "uMax = " << grid().uMax() << ", decades = " << decades << std::endl;
    // 1. Ensure that [decades(uM),decades(uP)] encloses [dM,dP],
    //    entirely or partially.
    if (decades<smartGrid.minEedfDecay)
    {
        // decades(uP) < dM: until decades(uP) >= dM, keep doubling uP,
        // set uMax=uP, solve and recalculate decades(uP)
        while (decades<smartGrid.minEedfDecay)
        {
            uP *= 2;
            if (grid().isUniform())
            {
                updateMaxEnergy(uP);
            } else
            {
                updateMaxEnergyNonuniform(uP);
            }
            solveSingle();
            decades = calcDecades(eedf[0],eedf[grid().nCells()-1]);
            std::cout << "uMax = " << grid().uMax() << ", decades = " << decades << std::endl;
        }
    }
    else if (decades>smartGrid.maxEedfDecay)
    {
        // decades(uM) > dP: until decades(uM) <= dP, keep dividing uM by 2,
        // set uMax=uM, solve and recalculate decades(uM)
        while (decades>smartGrid.maxEedfDecay)
        {
            uM /= 2;
            if (grid().isUniform())
            {
                updateMaxEnergy(uM);
            } else
            {
                updateMaxEnergyNonuniform(uM);
            }
            solveSingle();
            decades = calcDecades(eedf[0],eedf[grid().nCells()-1]);
            std::cout << "uMax = " << grid().uMax() << ", decades = " << decades << std::endl;
        }
    }
    else
    {
        // decades is already within the range: return early.
        std::cout << "Decades is within range: keeping the initial uMax." << std::endl;
        return;
    }
    // 2. If we reach this point, [decades(uM),decades(uP)] encloses
    //    [dM,dP], or part of it, and uMax is equal to uM or uP.
    //    while d(uMax) is not in [dM,dP], do a bisection of this interval:
    //      - set uMax = (um+uP)/2, solve and update decades
    //      - if d(uMax) is below dM, set uM=uMax
    //        else
    //        if d(uMax) is above dP, set uP=uMax
    while (decades<smartGrid.minEedfDecay || decades>smartGrid.maxEedfDecay)
    {
        // bisection step
        if (grid().isUniform())
        {
            updateMaxEnergy((uP+uM)/2);
        } else
        {
            updateMaxEnergyNonuniform((uP+uM)/2);
        }
        
        solveSingle();
        decades = calcDecades(eedf[0],eedf[grid().nCells()-1]);
        std::cout << "uMax = " << grid().uMax() << ", decades = " << decades << std::endl;
        if (decades<smartGrid.minEedfDecay)
        {
            uM = grid().uMax();
        }
        else
        {
            uP = grid().uMax();
        }
    }
    std::cout << "Final uMax = " << grid().uMax() << ", decades = " << decades << std::endl;
}

void ElectronKineticsBoltzmann::evaluateFirstAnisotropy()
{
    firstAnisotropy.setZero(grid().nCells());

    const double EoN = m_workingConditions->reducedFieldSI();
    const double WoN = m_workingConditions->reducedExcFreqSI();

    // 1. First fill firstAnisotropy with df/du.
    cellDerivative(firstAnisotropy,grid(),eedf);

    // make a copy. Below we add an ionization/attachment contribution
    Vector cellCrossSection(mixture.collision_data().totalCellCrossSection());

    if (ionizationOperator.includeNonConservativeIonization || attachmentOperator.includeNonConservativeAttachment)
    {
        if (growthModelType == GrowthModelType::temporal)
        {
            // The DC case. This implements Tejero2019, equation 3b:
            // f1 = -zeta(E/N)(1/Omega_PT)(df/du).

            // Calculate Omega_c. This is the first term that appears in
            // Tejero2019 equation 5a. It is defined in the text below 5b.
            cellCrossSection =
                cellCrossSection.array() + CIEff / (SI::gamma*grid().getCells().cwiseSqrt()).array();
            if (WoN == 0)
            {
                // DC case. zeta=1 and Omega_PT = Omega_c
                firstAnisotropy = -EoN * firstAnisotropy.cwiseQuotient(cellCrossSection);
            }
            else
            {
                // HF case. zeta=sqrt(2.) and Omega_PT includes the second term in 5a:
                // Omega_PT_i = Omega_c_i + (me/(2*e))(omega/N)^2 / (u*Omega_c_i)
                //   This can also be written as
                // Omega_PT_i = Omega_c_i + (omega/N)^2/(gamma^2*u*Omega_c_i)
                // (Note: that e*u is the energy in SI units.)
                // NOTE: EoN is the RMS field, but here we need the field amplitude. That explains the factor sqrt(2).
                firstAnisotropy = -EoN * std::sqrt(2.) * firstAnisotropy.array() /
                                  (cellCrossSection.array() +
                                   WoN * WoN / (SI::gamma*SI::gamma * grid().getCells().array() * cellCrossSection.array()));
            }
        }
        else if (growthModelType == GrowthModelType::spatial)
        {
            firstAnisotropy = -(alphaRedEff * eedf + EoN * firstAnisotropy).cwiseQuotient(cellCrossSection);
        }
    }
    else if (WoN == 0)
    {
        firstAnisotropy = -EoN * firstAnisotropy.cwiseQuotient(cellCrossSection);
    }
    else
    {
        // NOTE: EoN is the RMS field, but here we need the field amplitude. That explains the factor sqrt(2).
        firstAnisotropy =
            -EoN * std::sqrt(2.) * firstAnisotropy.array() /
            (cellCrossSection.array() + WoN * WoN / (SI::gamma*SI::gamma * grid().getCells().array() * cellCrossSection.array()));
    }
}

void ElectronKineticsBoltzmann::evaluatePower()
{
    const double EoN = m_workingConditions->reducedFieldSI();
    const double Tg = m_workingConditions->gasTemperature();
    // reset by an assignment of a default-constructed Power object
    power = Power();
    elasticOperator.evaluatePower(grid(),eedf,Tg,power.elasticNet,power.elasticGain,power.elasticLoss);
    if (carOperator.get())
    {
        carOperator->evaluatePower(grid(),eedf,Tg,power.carNet,power.carGain,power.carLoss);
    }
    fieldOperator.evaluatePower(grid(),eedf,EoN,power.field);
    // growth terms. Note that GrowthModelType::spatial also adds to power.field
    if (ionizationOperator.includeNonConservativeIonization || attachmentOperator.includeNonConservativeAttachment)
    {
        if (growthModelType == GrowthModelType::temporal)
        {
            power.eDensGrowth = -CIEff * getMeanEnergy(eedf,grid());
        }
        else if (growthModelType == GrowthModelType::spatial)
        {
            /** \todo field and correction are both of the form 'eedf*g'.
             *  Check that the following is correct. That requires that g_E and g_fieldSpatialGrowth
             *  have different dimensions (factor energy).
             */
            // Note that field is calculated by fieldOperator.evaluatePower, which already
            // adds the factor SI::gamma.
            const double correction = energyIntegral(grid(), interpolateNodalToCell(grid(), g_fieldSpatialGrowth), eedf);
            power.field -= SI::gamma * correction;

            const double powerMobility = - fNodegPrimeEnergyIntegral(grid(), grid().getNodes().array().pow(2) / mixture.collision_data().totalCrossSection().array(), eedf);
            const double powerDiffusion = energyIntegral(grid(), grid().getCells().cwiseAbs2().cwiseQuotient(mixture.collision_data().totalCellCrossSection()), eedf);
            power.eDensGrowth = - alphaRedEff * SI::gamma / 3. *
                                (powerMobility * m_workingConditions->reducedFieldSI() - powerDiffusion * alphaRedEff);
        }
    }

    if (eeOperator)
    {
        eeOperator->evaluatePower(grid(),eedf,power.electronElectronNet, power.electronElectronGain, power.electronElectronLoss);
    }
    // Evaluate power absorbed per electron at unit gas density due to in- and superelastic collisions.
    for (auto &cd : mixture.collision_data().data_per_gas())
    {
        cd.evaluatePower(ionizationOperator.ionizationOperatorType, eedf);
        power += cd.getPower();
    }
    /// \todo Change inelastic/superelastic with inelastic, use inelastic.forward, inelastic.backward.
    power.inelastic =
        power.excitation.forward
        + power.vibrational.forward
        + power.rotational.forward
        + power.ionization.forward
        + power.attachment.forward;
    power.superelastic =
        power.excitation.backward
        + power.vibrational.backward
        + power.rotational.backward;

    double totalGain = 0;
    /** \todo totalLoss is calculated, but never used.
     *  This produces a clang++ warning with -Wunused-but-set-variable.
     */
    double totalLoss = 0;

    /** \todo get rid of this; do 'for (double value : { ... this list ... })'
     *  in the line below.
     */
    double powerValues[14]{
        power.field,
        power.elasticGain, power.elasticLoss,
        power.carGain, power.carLoss,
        power.excitation.forward, power.excitation.backward,
        power.vibrational.forward, power.vibrational.backward,
        power.rotational.forward, power.rotational.backward,
        power.eDensGrowth,
        power.electronElectronGain, power.electronElectronLoss
    };

    for (double value : powerValues)
    {
        if (value > 0)
            totalGain += value;
        else
            totalLoss += value;
    }
    power.balance =
        power.field
        + power.elasticNet
        + power.carNet
        + power.inelastic
        + power.superelastic
        + power.eDensGrowth
        + power.electronElectronNet;
    power.relativeBalance = std::abs(power.balance) / totalGain;
    power.reference = totalGain;
}

void ElectronKineticsBoltzmann::evaluateSwarmParameters()
{
    const bool nonConservative = ionizationOperator.includeNonConservativeIonization
                              || attachmentOperator.includeNonConservativeAttachment;

    /* First we calculate the DC mobilities and diffusion coefficients. These
     * are based on Omega_SST (spatial growth) or Omega_c (temporal growth)
     * (manual v2.2 eqns. 46a,b). From eqns. 11a,b) we learn that Omega_SST is
     * just sigma_c; for temporal growth we calculate Omega_c by adding the
     * ionization term.
     */
    /** \todo For spatial growth (or no growth), the copy is not needed; we
     *  could use totalCellCrossSection() directly.
     */
    Vector cellCS(mixture.collision_data().totalCellCrossSection());

    if (growthModelType == GrowthModelType::temporal && nonConservative)
    {
        cellCS.array() += (CIEff/SI::gamma) / grid().getCells().cwiseSqrt().array();
    }
    const Vector D0 = grid().getCells().array() / (3. * cellCS).array();
    const Vector D0Nodes = grid().getNodes().array() / (3. * mixture.collision_data().totalCrossSection()).array();

    swarmParameters.redDiffCoeff = SI::gamma*energyIntegral(grid(),D0,eedf);
    swarmParameters.redDiffCoeffEnergy = SI::gamma*energyIntegral(grid(),grid().getCells().cwiseProduct(D0),eedf);
    swarmParameters.redMobCoeff = -SI::gamma*fNodegPrimeEnergyIntegral(grid(),D0Nodes,eedf);
    swarmParameters.redMobilityEnergy = -SI::gamma*fNodegPrimeEnergyIntegral(grid(),grid().getNodes().cwiseProduct(D0Nodes),eedf);

    /* If WoN !=0, we also calculate the HF mobility (real and imaginary parts).
     * Note that WoN implies the temporal growth case, so cellCS represents
     * Omega_c (we just added the ionization term). First we construct Omega_PT
     * by copying Omega_c a and adding the WoN term (eq. 11a). Next we calculate
     * Re{mu_HF} and Im{mu_HF} from Omega_PT and Omega_c using eqns. 49a,b).
     */
    const double WoN = m_workingConditions->reducedExcFreqSI();
    if (WoN>0.0)
    {
        const Vector& OmegaC = cellCS;
        const Vector OmegaPT = OmegaC + ( WoN * WoN / (SI::gamma*SI::gamma)) * (grid().getCells().cwiseProduct(OmegaC)).cwiseInverse();
        double muHFRe = -SI::gamma / 3. * fgPrimeEnergyIntegral(grid(),grid().getCells().cwiseQuotient(OmegaPT),eedf);
        double muHFIm = (+SI::gamma / 3. * WoN/SI::gamma) *
                      fgPrimeEnergyIntegral(grid(),grid().getCells().cwiseSqrt().cwiseQuotient(OmegaPT).cwiseQuotient(OmegaC),eedf);
        swarmParameters.redMobilityHF = { muHFRe, muHFIm };
    }

    if (growthModelType == GrowthModelType::spatial && nonConservative)
    {
        swarmParameters.driftVelocity = -swarmParameters.redDiffCoeff * alphaRedEff +
                                        swarmParameters.redMobCoeff * m_workingConditions->reducedFieldSI();
    }
    else
    {
        /// \todo Does this make sense in case of a HF field?
        swarmParameters.driftVelocity = swarmParameters.redMobCoeff * m_workingConditions->reducedFieldSI();
    }

    double totalIonRateCoeff = 0., totalAttRateCoeff = 0.;

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (const auto &collision : cd.collisions(CollisionType::ionization))
        {
            totalIonRateCoeff += collision->getTarget()->delta() * collision->ineRateCoeff();
        }

        for (const auto &collision : cd.collisions(CollisionType::attachment))
        {
            totalAttRateCoeff += collision->getTarget()->delta() * collision->ineRateCoeff();
        }
    }
    swarmParameters.redTownsendCoeff = totalIonRateCoeff / swarmParameters.driftVelocity;
    swarmParameters.redAttCoeff = totalAttRateCoeff / swarmParameters.driftVelocity;

    swarmParameters.meanEnergy = getMeanEnergy(eedf,grid());

    if (swarmParameters.meanEnergy < 0)
    {
        Log<Message>::Warning("Negative mean electron energy ", swarmParameters.meanEnergy,
                              "eV encountered at E/N = ", m_workingConditions->reducedField(),
                              "Td. Make sure to use an appropriate grid resolution and maximum energy.");
    }

    swarmParameters.characEnergy = swarmParameters.redDiffCoeff / swarmParameters.redMobCoeff;

    swarmParameters.Te = 2. / 3. * swarmParameters.meanEnergy;

    // TODO: is this correct? (simulations after the first will have a different value for Te).
    m_workingConditions->updateElectronTemperature(swarmParameters.Te);
}



ElectronKineticsPrescribed::ElectronKineticsPrescribed(const std::filesystem::path &basePath,const json_type &cnf, WorkingConditions *workingConditions)
: ElectronKinetics(basePath,cnf,workingConditions),
  shapeParameter(cnf.at("shapeParameter").get<double>())
{
    if (getIonizationOperatorType(cnf.at("ionizationOperatorType")) != IonizationOperatorType::conservative)
    {
        throw std::runtime_error("ionizationOperatorType must be 'conservative' for EEDF type 'Prescribed'.");
    }
    grid().updatedMaxEnergy.addListener(&ElectronKineticsPrescribed::evaluateMatrix, this);
    m_workingConditions->updatedReducedField.addListener(&ElectronKineticsPrescribed::evaluateFieldOperator, this);
    /// \todo Do this here? Not all parameters may have been set at this point.
    evaluateMatrix();
}

void ElectronKineticsPrescribed::evaluateFieldOperator()
{
    const double WoN = m_workingConditions->reducedExcFreqSI();
    const double dummyCIEff = 0.0;
    fieldOperator.evaluate(grid(),mixture.collision_data().totalCrossSection(),WoN,dummyCIEff);
}

void ElectronKineticsPrescribed::evaluateMatrix()
{
    mixture.collision_data().evaluateTotalAndElasticCS(grid());

    elasticOperator.evaluate(grid(),mixture.collision_data().elasticCrossSection());

    /** This is most probably not correct. For the Prescribed EEDF mode, an equivalent
     *  field must be determined from the enery balance, whereas evaluateFieldOperator()
     *  uses the prescribed E/N.
     */
    evaluateFieldOperator();

    if (carOperator.get())
    {
        carOperator->evaluate(grid());
    }

    inelasticOperator.evaluateInelasticOperators(grid(),mixture);
}

void ElectronKineticsPrescribed::doSolve()
{
    /** This implements the logic from the MATLAB code in Prescribed.m:
     *  1. Inspect smartGrid(). When active, calculate the number of decades as the
     *     mean value of the minumum and maximum numbers of dcades and calculate
     *     u_max from the result. Update the grid() accordingly.
     *  2. Assign values to eedf using the shape parameter and the present value of Te
     *     that is obtained from the WorkingConditions.
     *  3. Calculate relevant derived values (power, rate coefficients, swarm parameters).
     */
    const double Te=m_workingConditions->electronTemperature();
    if (grid().smartGrid())
    {
        const double g = shapeParameter;
        const double gamma_3_2g = std::tgamma(3/(2*g));
        const double gamma_5_2g = std::tgamma(5/(2*g));
        double decades = 0.5*(grid().smartGrid()->minEedfDecay + grid().smartGrid()->maxEedfDecay);
        double maxEnergy = std::pow(decades/std::log10(std::exp(1)),1/g)*1.5*Te*gamma_3_2g/gamma_5_2g;
        updateMaxEnergy(maxEnergy);
    }
    // evaluate the EEDF
    makePrescribedEDF(eedf,grid(),shapeParameter,Te);

    evaluatePower();

    // evaluate derived data
    mixture.collision_data().evaluateRateCoefficients(grid(),eedf);
    evaluateSwarmParameters();

    // nullptr: firstAnisotropy is not available
    obtainedNewEedf.emit(grid(), eedf, *m_workingConditions, power, mixture.collision_data(), swarmParameters,
                         nullptr);
}


void ElectronKineticsPrescribed::evaluatePower()
{
    const double Tg = m_workingConditions->gasTemperature();

    // reset by an assignment of a default-constructed Power object
    power = Power();

    elasticOperator.evaluatePower(grid(),eedf,Tg,power.elasticNet,power.elasticGain,power.elasticLoss);

    if (carOperator.get())
    {
        carOperator->evaluatePower(grid(),eedf,Tg,power.carNet,power.carGain,power.carLoss);
    }

    // Evaluate power absorbed per electron at unit gas density due to in- and superelastic collisions.
    for (auto &cd : mixture.collision_data().data_per_gas())
    {
        // note: in the constructor it is tested that this was specified by the user.
        // other types do not make sense for the prescribed eedf case, in which also
        // no growth model is assumed.
        const IonizationOperatorType ionizationOperatorType = IonizationOperatorType::conservative;
        cd.evaluatePower(ionizationOperatorType, eedf);
        power += cd.getPower();
    }

    /// \todo Change inelastic/superelastic with inelastic, use inelastic.forward, inelastic.backward.
    power.inelastic =
        power.excitation.forward
        + power.vibrational.forward
        + power.rotational.forward
        + power.ionization.forward
        + power.attachment.forward;
    power.superelastic =
        power.excitation.backward
        + power.vibrational.backward
        + power.rotational.backward;

    power.balance =
        /* skip power.field, see below... */
        + power.elasticNet
        + power.carNet
        + power.inelastic
        + power.superelastic
        + power.eDensGrowth
        + power.electronElectronNet;

    /* next, calculate E/N from the requirement that power.balance==0.
     * we have P_tot = P_field + P_other, and P_field(E/N) = C*(E/N)^2 for some C.
     * 1) This C can be calculated as P_field(E/N=1). We then find that
     * 2) P_tot = 0 => P_field = -P_other, and
     * 3) P_field = P_field(E/N=1)*(E/N)^2 => E/N = sqrt(P_field/P_field(E/N=1)).
     */
    // 1. Calculate P_1 := P_field(E/N=1)
    const double EoN = m_workingConditions->reducedFieldSI();
    const double WoN = m_workingConditions->reducedExcFreqSI();
    const double dummyCIEff= 0.0;
    fieldOperator.evaluate(grid(),mixture.collision_data().totalCrossSection(),WoN,dummyCIEff);
    double P_1;
    fieldOperator.evaluatePower(grid(),eedf,EoN,P_1);
    // 2. we calculate P_E such that the power balance will be satisfied
    power.field = -power.balance;
    // 3. E_N = sqrt(P_field/P_field(E/N=1)).
    //    set this E/N value in the workingCOndition object, so it will be printed correctly.
    m_workingConditions->updateReducedField(std::sqrt(power.field/P_1)/SI::Townsend);

    // by construction, the power balance will now be satisfied.
    power.balance = 0.0;
    power.relativeBalance = 0.0;
    power.reference = 0.0;
}

void ElectronKineticsPrescribed::evaluateSwarmParameters()
{
    const Vector& cellCS(mixture.collision_data().totalCellCrossSection());
    const Vector D0 = grid().getCells().array() / (3. * cellCS).array();

    swarmParameters.redDiffCoeff = SI::gamma*energyIntegral(grid(),D0,eedf);
    swarmParameters.redDiffCoeffEnergy = SI::gamma*energyIntegral(grid(),grid().getCells().cwiseProduct(D0),eedf);
    swarmParameters.redMobCoeff = -SI::gamma*fgPrimeEnergyIntegral(grid(),D0,eedf);
    swarmParameters.redMobilityEnergy = -SI::gamma*fgPrimeEnergyIntegral(grid(),grid().getCells().cwiseProduct(D0),eedf);

    const double WoN = m_workingConditions->reducedExcFreqSI();
    if (WoN>0.0)
    {
        /* NOTE: in the prescribed case, only conservative processes should be
         * taken into account. In that case CIEff==0 and OmegaC==cellCS. Otherwise,
         * the code below is equal to the code in the function
         * ElectronKineticsBoltzmann::evaluateSwarmParameters. We highlight that
         * fact by using OmegaC, but making it just a reference to cellCS.
         */
        const Vector& OmegaC = cellCS;
        const Vector OmegaPT = OmegaC + ( WoN * WoN / (SI::gamma*SI::gamma)) * (grid().getCells().cwiseProduct(OmegaC)).cwiseInverse();
        double muHFRe = -SI::gamma / 3. * fgPrimeEnergyIntegral(grid(),grid().getCells().cwiseQuotient(OmegaPT),eedf);
        double muHFIm = (+SI::gamma / 3. * WoN/SI::gamma) *
                      fgPrimeEnergyIntegral(grid(),grid().getCells().cwiseSqrt().cwiseQuotient(OmegaPT).cwiseQuotient(OmegaC),eedf);
        swarmParameters.redMobilityHF = { muHFRe, muHFIm };
    }
    swarmParameters.driftVelocity = swarmParameters.redMobCoeff * m_workingConditions->reducedFieldSI();

    double totalIonRateCoeff = 0., totalAttRateCoeff = 0.;

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (const auto &collision : cd.collisions(CollisionType::ionization))
        {
            totalIonRateCoeff += collision->getTarget()->delta() * collision->ineRateCoeff();
        }

        for (const auto &collision : cd.collisions(CollisionType::attachment))
        {
            totalAttRateCoeff += collision->getTarget()->delta() * collision->ineRateCoeff();
        }
    }

    swarmParameters.redTownsendCoeff = totalIonRateCoeff / swarmParameters.driftVelocity;
    swarmParameters.redAttCoeff = totalAttRateCoeff / swarmParameters.driftVelocity;

    swarmParameters.meanEnergy = getMeanEnergy(eedf,grid());

    swarmParameters.characEnergy = swarmParameters.redDiffCoeff / swarmParameters.redMobCoeff;

    /* Just assign the working condition Te (which is input), as is also done in MATLAB.
     * NOTE that in v1.0.0 of the MATLAB code (Code/PrescribedEedf.m line 365) there is
     * a comment that suggests the calculation Te = 2./3.*swarmParameters.meanEnergy. This
     * comment disappeared in version 2.2.0. Also in the C++ version the calculation was
     * implemented as
     *   swarmParameters.Te = 2. / 3. * swarmParameters.meanEnergy;
     * until 9f0200138acf2ca3006fb0a27af524002ce41d41 (including)
     * In an ideal world, these are the same, but there will be differences because of the
     * discretization of the EEDF and the truncation to some u_max. It would be interesting
     * to print both values, the discrepancy is a measure for the discretization and energy
     * truncation errors.
     */
    swarmParameters.Te = m_workingConditions->electronTemperature();
}

} // namespace loki
