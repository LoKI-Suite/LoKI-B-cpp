//
// Created by daan on 13-5-19.
//

#include "LoKI-B/ElectronKinetics.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/EedfUtilities.h"
#include <chrono>
#include <cmath>

namespace fs = std::filesystem;

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

ElectronKinetics::ElectronKinetics(const fs::path &basePath, const ElectronKineticsSetup &setup, WorkingConditions *workingConditions)
    : m_workingConditions(workingConditions),
    m_grid(setup.numerics.energyGrid),
    mixture(basePath, &grid(), setup, workingConditions),
    fieldOperator(grid()),
    inelasticOperator(grid()),
    eedf(grid().nCells())
{
    initialize();
}

ElectronKinetics::ElectronKinetics(const fs::path &basePath, const json_type &cnf, WorkingConditions *workingConditions)
    : m_workingConditions(workingConditions),
    m_grid(Grid::fromConfig(cnf.at("numerics").at("energyGrid"))),
    mixture(basePath, &grid(), cnf, workingConditions),
    fieldOperator(grid()),
    inelasticOperator(grid()),
    eedf(grid().nCells())
{
    initialize();
}

void ElectronKinetics::initialize()
{
    if (!mixture.CARGases().empty())
    {
        carOperator.reset(new CAROperator(mixture.CARGases()));
    }
}

void ElectronKinetics::updateMaxEnergy(double uMax)
{
    m_grid.updateMaxEnergy(uMax);
}

void ElectronKinetics::solve()
{
    doSolve();
}

ElectronKineticsBoltzmann::ElectronKineticsBoltzmann(const std::filesystem::path &basePath, const ElectronKineticsSetup &setup, WorkingConditions *workingConditions)
: ElectronKinetics(basePath,setup,workingConditions),
    ionizationOperator(setup.ionizationOperatorType),
    eeOperator(setup.includeEECollisions ? new ElectronElectronOperator(grid()) : nullptr),
    fieldMatrixSpatGrowth(grid().nCells(), grid().nCells()),
    ionSpatialGrowthD(grid().nCells(), grid().nCells()),
    ionSpatialGrowthU(grid().nCells(), grid().nCells()),
    fieldMatrixTempGrowth(grid().nCells(), grid().nCells()),
    ionTemporalGrowth(grid().nCells(), grid().nCells())
{
    this->mixingParameter = setup.numerics.nonLinearRoutines.mixingParameter;
    this->maxEedfRelError = setup.numerics.nonLinearRoutines.maxEedfRelError;
    this->maxPowerBalanceRelError = setup.numerics.maxPowerBalanceRelError;
    this->growthModelType = setup.growthModelType;
    initialize();
}

ElectronKineticsBoltzmann::ElectronKineticsBoltzmann(const std::filesystem::path &basePath, const json_type &cnf, WorkingConditions *workingConditions)
: ElectronKinetics(basePath, cnf,workingConditions),
    ionizationOperator(getIonizationOperatorType(cnf.at("ionizationOperatorType"))),
    eeOperator(cnf.at("includeEECollisions") ? new ElectronElectronOperator(grid()) : nullptr),
    fieldMatrixSpatGrowth(grid().nCells(), grid().nCells()),
    ionSpatialGrowthD(grid().nCells(), grid().nCells()),
    ionSpatialGrowthU(grid().nCells(), grid().nCells()),
    fieldMatrixTempGrowth(grid().nCells(), grid().nCells()),
    ionTemporalGrowth(grid().nCells(), grid().nCells())
{
    this->mixingParameter = cnf.at("numerics").at("nonLinearRoutines").at("mixingParameter");
    this->maxEedfRelError = cnf.at("numerics").at("nonLinearRoutines").at("maxEedfRelError");
    this->maxPowerBalanceRelError = cnf.at("numerics").at("maxPowerBalanceRelError");
    this->growthModelType = getGrowthModelType(cnf.at("growthModelType"));
    initialize();
}

void ElectronKineticsBoltzmann::initialize()
{
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
    const double EoN = m_workingConditions->reducedFieldSI();
    const double WoN = m_workingConditions->reducedExcFreqSI();
    fieldOperator.evaluate(grid(),mixture.collision_data().totalCrossSection(),EoN,WoN,fieldMatrix);
}

void ElectronKineticsBoltzmann::evaluateMatrix()
{
    const double Tg = m_workingConditions->gasTemperature();
    mixture.collision_data().evaluateTotalAndElasticCS(grid());

    elasticOperator.evaluate(grid(),mixture.collision_data().elasticCrossSection(),Tg,elasticMatrix);

    evaluateFieldOperator();

    if (carOperator.get())
    {
        carOperator->evaluate(grid(),Tg,CARMatrix);
    }

    inelasticOperator.evaluateInelasticOperators(grid(),mixture);

    if (mixture.collision_data().hasCollisions(CollisionType::ionization))
        ionizationOperator.evaluateIonizationOperator(grid(),mixture);

    if (mixture.collision_data().hasCollisions(CollisionType::attachment))
        attachmentOperator.evaluateAttachmentOperator(grid(),mixture);

    /** \todo see the comments about superElasticThresholds in the header file.
    // Sort and erase duplicates.
    std::sort(superElasticThresholds.begin(), superElasticThresholds.end());
    superElasticThresholds.erase(unique(superElasticThresholds.begin(), superElasticThresholds.end()),
                                 superElasticThresholds.end());
    */
    if (eeOperator)
    {
        eeOperator->initialize(grid());
    }
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

void ElectronKineticsBoltzmann::invertLinearMatrix()
{
    // Here the Conservative ionization and attachment matrices are added.
    // Otherwise, the sum is the same as in solveSpatialGrowthMatrix
    // (note that in solveTemporalGrowthMatrix fieldMatrix is not added).
    boltzmannMatrix
        = elasticMatrix
        + fieldMatrix
        + inelasticOperator.inelasticMatrix
        + ionizationOperator.ionConservativeMatrix
        + attachmentOperator.attachmentConservativeMatrix;
    if (carOperator)
    {
        boltzmannMatrix += CARMatrix;
    }

    invertMatrix(boltzmannMatrix);
}

void ElectronKineticsBoltzmann::invertMatrix(Matrix &matrix)
{
//#define LOKIB_TIME_INVERT_MATRIX
#ifdef LOKIB_TIME_INVERT_MATRIX
    auto begin = std::chrono::high_resolution_clock::now();
#endif

    if (!inelasticOperator.hasSuperelastics)
    {
        eedf.setZero();
        eedf[0] = 1.;

        if (grid().isUniform())
        {
            matrix.row(0) = grid().getCells().cwiseSqrt() * grid().du();
        } else
        {
            matrix.row(0) = grid().getCells().cwiseSqrt().cwiseProduct(grid().duCells());
        }
        /** We can also just make eedf[0]=1 by implementing the above and:
         *  matrix.row(0).setZero(); matrix(0,0)=1.0;
         *  The advantage is that the first row only has the diagonal,
         *  over-all the matrix has a better sparsity pattern.
         *  Then we can do the normalization afterwards.
         */

        LinAlg::hessenberg(matrix.data(), eedf.data(), grid().nCells());
    }
    else
    {
        /** \todo Explain the idea that is outlined below. Is it correct that p is not used
         *        after the call to LinAlg::hessenbergReductionPartialPiv? Is this icomplete?
         */
        // TODO: Find a way to distinguish when to use LU and when to use Hessenberg reduction.

        // HESSENBERG WITH PARTIAL PIVOTING
        //            auto *p = new uint32_t[grid().nCells()];
        //
        //            for (Grid::Index i = 0; i < grid().nCells(); ++i)
        //                p[i] = i;
        //
        //            LinAlg::hessenbergReductionPartialPiv(matrix.data(), &superElasticThresholds[0], p,
        //              grid().nCells(), superElasticThresholds.size());
        //
        //            eedf.setZero();
        //            eedf[0] = 1.;
        //
        //            matrix.row(0) = grid().getCells().cwiseSqrt() * grid().du();
        //
        //            LinAlg::hessenberg(matrix.data(), eedf.data(), grid().nCells());
        //
        //            delete[] p;

        // LU DECOMPOSITION
        Vector b = Vector::Zero(grid().nCells());
        b[0] = 1;

        if (grid().isUniform())
        {
            matrix.row(0) = grid().getCells().cwiseSqrt() * grid().du();
        } else
        {
            matrix.row(0) = grid().getCells().cwiseSqrt().cwiseProduct(grid().duCells());
        }
        eedf = matrix.partialPivLu().solve(b);
    }

    /** \todo It seems that the normaization is superfluous, since the normalization condition
     *        is already part of the system (first row of A, first element of b). One could
     *        decide to change the first equation into eedf[0] = 1 and do the normalization
     *        afterwards. That prevents a fully populated first row of the system matrix
     *        (better sparsity pattern).
     */
    // std::cout << "NORM: " <<  eedf.dot(grid().getCells().cwiseSqrt() * grid().du()) << std::endl;
    if (grid().isUniform())
        {
            eedf /= eedf.dot(grid().getCells().cwiseSqrt() * grid().du());
        } else
        {
            eedf /= eedf.dot(grid().getCells().cwiseSqrt().cwiseProduct(grid().duCells()));
        }

#ifdef LOKIB_TIME_INVERT_MATRIX
    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "Inverted matrix elapsed time = "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;
#endif
}

void ElectronKineticsBoltzmann::solveSpatialGrowthMatrix()
{
    const double EoN = m_workingConditions->reducedFieldSI();

    Vector cellTotalCrossSection(grid().nCells());

    for (Grid::Index i = 0; i < grid().nCells(); ++i)
    {
        cellTotalCrossSection[i] = .5 * (mixture.collision_data().totalCrossSection()[i] + mixture.collision_data().totalCrossSection()[i + 1]);
    }
    boltzmannMatrix
        = elasticMatrix
        + fieldMatrix
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

    /** \todo The name is incorrect. This is not the integrand since it already includes
     *        du and does not yet have the factor f(u).
     */
    Vector integrandCI;
    if (grid().isUniform())
        {
            integrandCI = SI::gamma*grid().du() * Vector::Ones(grid().nCells()).transpose() *
                         (ionizationOperator.ionizationMatrix + attachmentOperator.attachmentMatrix);
        } else
        {
            integrandCI = SI::gamma*grid().duCells().transpose() *
                         (ionizationOperator.ionizationMatrix + attachmentOperator.attachmentMatrix);
        }

    double CIEffNew = eedf.dot(integrandCI);
    double CIEffOld = CIEffNew / 3;
    /** \todo Where does the division by 3 come from? Without that, CIEff
     *  is calculated below as a weighted average of old and new values
     *  (underrelaxation). I cannot interpret the equation with the additional /3.
     *  NOTE that this is just initialization; in the while(!converged) loop
     *  we do not have such factor 1/3.
     *  Why do we need an 'old' value on entry of that loop?
     */
    CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

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
     * Subsequently, we evaluate D_eN=+gamma*int_0^infty D0(u)f(u)du (eqn. 19a),
     * which is straightforward.
     *
     * For the mobility we have mu_eE=-gamma*int_0^infty D0(u)[df/du]du (19b).
     * which is handled quite differently. It appears as if in the code below
     * a partial integration is carried out: since on the boundaries u=0 and
     * u=infty we have uf(u)=0 this gives int D(u)f'(u)du = -int D'(u)f(u)du.
     * In the code below we define U0=-D'(u), which results in the expression
     * mu_eN=-gamma*int_0^infty U0(u)f(u) du. This is what you find below.
     */
    /** \todo The calculation of U0 on the boundaries needs to be explained.
     *  These seem to be incorrect. In the approximation D'[i]=(D[i+1]-D[i-1])/(2*du)
     *  the problematic terms i-1 (at i=0) and i+i (at i=Nc) are simply omitted.
     *  That results in uncontrolled discretization errors (in practice small,
     *  if you have enough grid() points).
     */
    /** \todo cellTotalCrossSection is the result of interpolation. See if
     *  it may be better to use the original CS at the nodes in parts of these
     *  calculations. In particular, the calculation of U will be much more
     *  straighforward if we do a node->cell interpolation, since also the
     *  first and last cells have two node-neighbours. Does does also
     *  influence the accuracy with which invariants are reproduced (for
     *  example: the Einstein relation or the characteristic temperature in
     *  case of a Maxwellian eedf)? Check such things first.
     */
    const Vector D0 = grid().getCells().array() / (3. * cellTotalCrossSection).array();
    /** \todo Document/explain which equation is discretized here.
     */
    Vector U0sup(grid().nCells());
    Vector U0inf(grid().nCells());
    U0sup[0] = 0.;
    U0inf[grid().nCells() - 1] = 0.;
    for (Grid::Index j = 0; j < grid().nCells(); ++j)
    {
        if (j != 0)
            U0sup[j] = EoN / (2. * grid().du()) * D0[j - 1];

        if (j != grid().nCells() - 1)
            U0inf[j] = -EoN / (2. * grid().du()) * D0[j + 1];
    }
    const Vector U0 = U0sup + U0inf;

    // This is 33a from \cite Manual_1_0_0
    double ND  =   SI::gamma * grid().du() * D0.dot(eedf);
    /* This is 33b from \cite Manual_1_0_0, multiplied with E/N.
     * Note that the factor E/N is part of U0sup, U0inf.
     */
    double muE = - SI::gamma * grid().du() * U0.dot(eedf);

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
    /** The following line implements the solution of equation 22 of
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

        for (Grid::Index k = 0; k < grid().nCells(); ++k)
        {
            // note: fieldMatrixSpatGrowth represents: (alphaEffNew/N)*(E/N)*d(D^0*f0)/du
            fieldMatrixSpatGrowth.coeffRef(k, k) = (g_fieldSpatialGrowth[k + 1] - g_fieldSpatialGrowth[k]) / (2*grid().du());
            // note: this is (alphaEffNew/N)^2*D0[k]
            ionSpatialGrowthD.coeffRef(k, k) = alphaRedEffNew * alphaRedEffNew * D0[k];
            boltzmannMatrix(k, k) = baseDiag[k] + fieldMatrixSpatGrowth.coeff(k, k) +
                                                           ionSpatialGrowthD.coeff(k, k);

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
#define USE_D0_FOR_ionSpatialGrowthU 1
            if (k > 0)
            {
                fieldMatrixSpatGrowth.coeffRef(k, k - 1) = -g_fieldSpatialGrowth[k] / (2*grid().du());
#if USE_D0_FOR_ionSpatialGrowthU
                ionSpatialGrowthU.coeffRef(k, k - 1) = -alphaRedEffNew*EoN*D0[k] / (2.*grid().du());
#else
                ionSpatialGrowthU.coeffRef(k, k - 1) = alphaRedEffNew * U0inf[k - 1];
#endif
                boltzmannMatrix(k, k - 1) = baseSubDiag[k] + fieldMatrixSpatGrowth.coeff(k, k - 1) +
                                                                      ionSpatialGrowthU.coeff(k, k - 1);
            }

            if (k < grid().nCells() - 1)
            {
                fieldMatrixSpatGrowth.coeffRef(k, k + 1) = g_fieldSpatialGrowth[k + 1] / (2*grid().du());
#if USE_D0_FOR_ionSpatialGrowthU
                ionSpatialGrowthU.coeffRef(k, k + 1) = +alphaRedEffNew*EoN*D0[k] / (2.*grid().du());
#else
                ionSpatialGrowthU.coeffRef(k, k + 1) = alphaRedEffNew * U0sup[k + 1];
#endif
                boltzmannMatrix(k, k + 1) = baseSupDiag[k] + fieldMatrixSpatGrowth.coeff(k, k + 1) +
                                                                      ionSpatialGrowthU.coeff(k, k + 1);
            }
        }

        Vector eedfNew = eedf;

        invertMatrix(boltzmannMatrix);

        CIEffOld = CIEffNew;
        CIEffNew = eedf.dot(integrandCI);

        CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

        ND  =   SI::gamma * grid().du() * D0.dot(eedf);
        muE = - SI::gamma * grid().du() * U0.dot(eedf);

        /** \todo See the notes just above this loop for a note about the discontinuity.
         */
        const double discriminant_new = muE * muE - 4 * CIEffNew * ND;

        alphaRedEffOld = alphaRedEffNew;
        alphaRedEffNew = (discriminant_new < 0) ? CIEffNew / muE : (muE - std::sqrt(discriminant_new)) / (2 * ND);

        alphaRedEffNew = mixingParameter * alphaRedEffNew + (1 - mixingParameter) * alphaRedEffOld;

        if (((alphaRedEffNew == 0 || std::abs(alphaRedEffNew - alphaRedEffOld) / alphaRedEffOld < 1.e-10) &&
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
        ++iter;
    }

    std::cerr << "Spatial growth routine converged in: " << iter << " iterations.\n";

    alphaRedEff = alphaRedEffOld;
    CIEff = CIEffOld;
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
    const Vector integrandCI = (SI::gamma * grid().du())
                    * Vector::Ones(grid().nCells()).transpose()
                    * (ionizationOperator.ionizationMatrix + attachmentOperator.attachmentMatrix);
    double CIEffNew = eedf.dot(integrandCI);
    double CIEffOld = CIEffNew / 3.;
    CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

    Vector eedfNew(grid().nCells());

    bool hasConverged = false;
    uint32_t iter = 0;

    while (!hasConverged)
    {
 //       Log<Message>::Notify("Iteration ", iter);

        // CIEff is <nu_eff>/N, so growthFactor = <nu_eff>/(N*gamma)
        const long double growthFactor = CIEffNew / SI::gamma;

        g_fieldTemporalGrowth.resize(grid().getNodes().size());
        g_fieldTemporalGrowth[0] = 0.;
        for (Grid::Index i=1; i!= g_fieldTemporalGrowth.size()-1; ++i)
        {
            // totalCSI = Omega_PT, See \cite Tejero2019 5a or the Manual 2.2.0, below eq. 11b:
            const double totalCSI = mixture.collision_data().totalCrossSection()[i] + growthFactor / std::sqrt(grid().getNode(i));
            /* The following g corresponds to the G_E in \cite Tejero2019 equation 6a,
             * with f^1(u) as in equation 3b
             * In the Manual 2.2.0: G_E as in 12a, f^1(u) as in 7b
             */
            const double OmegaPT = totalCSI + ( WoN * WoN / (SI::gamma*SI::gamma)) / (grid().getNode(i)*totalCSI);
            g_fieldTemporalGrowth[i] = (EoN * EoN / 3) * grid().getNode(i) / OmegaPT;
        }
        g_fieldTemporalGrowth[g_fieldTemporalGrowth.size() - 1] = 0.;

        const double sqrStep = grid().du() * grid().du();

        for (Grid::Index k = 0; k < grid().nCells(); ++k)
        {
            fieldMatrixTempGrowth.coeffRef(k, k) = -(g_fieldTemporalGrowth[k] + g_fieldTemporalGrowth[k + 1]) / sqrStep;
            // Manual 2.2.0, 7a (with an minus sign because all terms are negated):
            ionTemporalGrowth.coeffRef(k, k) = -growthFactor * std::sqrt(grid().getCell(k));
            boltzmannMatrix(k, k) = baseDiag[k] + fieldMatrixTempGrowth.coeff(k, k) +
                                                           ionTemporalGrowth.coeff(k, k);

            if (k > 0)
            {
                fieldMatrixTempGrowth.coeffRef(k, k - 1) = g_fieldTemporalGrowth[k] / sqrStep;
                boltzmannMatrix(k, k - 1) = baseSubDiag[k] + fieldMatrixTempGrowth.coeff(k, k - 1);
            }

            if (k < grid().nCells() - 1)
            {
                fieldMatrixTempGrowth.coeffRef(k, k + 1) = g_fieldTemporalGrowth[k + 1] / sqrStep;
                boltzmannMatrix(k, k + 1) = baseSupDiag[k] + fieldMatrixTempGrowth.coeff(k, k + 1);
            }

        }

        eedfNew = eedf;

        invertMatrix(boltzmannMatrix);

        CIEffOld = CIEffNew;
        CIEffNew = eedf.dot(integrandCI);
        CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

        if (((CIEffNew == 0 || std::abs(CIEffNew - CIEffOld) / CIEffOld < 10e-10) &&
             maxRelDiff(eedfNew,eedf) < maxEedfRelError) ||
            iter > 150)
        {
            hasConverged = true;

            if (iter > 150 && !eeOperator)
                Log<Message>::Warning("Iterative temporal growth scheme did not converge.");
        }

        ++iter;
    }

 //   std::cerr << "Temporal growth routine converged in: " << iter << " iterations.\n";

    CIEff = CIEffOld;
}

void ElectronKineticsBoltzmann::solveEEColl()
{
    assert(eeOperator);
    // Splitting all possible options for the best performance.

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
                + fieldMatrix
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
                + fieldMatrixTempGrowth;
        }
    }
    else
    {
        boltzmannMatrix
                = ionizationOperator.ionConservativeMatrix
                + attachmentOperator.attachmentConservativeMatrix
                + elasticMatrix
                + inelasticOperator.inelasticMatrix
                + fieldMatrix;
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

        const double ratio = std::abs(power.electronElectron / power.reference);

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
    invertLinearMatrix();

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
        updateMaxEnergy(grid().uMax() * (1 + smartGrid.updateFactor));
        solveSingle();
        decades = calcDecades(eedf[0],eedf[grid().nCells()-1]);
        //std::cout << "decades: " << decades << ", uMax: " << grid().uMax() << std::endl;
    }

    while (decades > smartGrid.maxEedfDecay)
    {
        updateMaxEnergy(grid().uMax() / (1 + smartGrid.updateFactor));
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
            updateMaxEnergy(uP);
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
            updateMaxEnergy(uM);
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
        updateMaxEnergy((uP+uM)/2);
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
    const Grid::Index n = grid().nCells();

    // 1. First fill firstAnisotropy with df/du.
    if (grid().isUniform())
    {
        firstAnisotropy[0] = (eedf[1] - eedf[0]) / grid().du();
        firstAnisotropy[n - 1] = (eedf[n - 1] - eedf[n - 2]) / grid().du();
        firstAnisotropy.segment(1, n - 2) = (eedf.segment(2, n - 2) - eedf.segment(0, n - 2)) / (2 * grid().du());
    } else
    {
        firstAnisotropy[0] = (eedf[1] - eedf[0]) / grid().duNode(1);
        firstAnisotropy[n - 1] = (eedf[n - 1] - eedf[n - 2]) / grid().duNode(n-1);
        firstAnisotropy.segment(1, n - 2) = (eedf.segment(2, n - 2) - eedf.segment(0, n - 2)).cwiseQuotient(grid().duNodes().segment(1, n - 2) + grid().duNodes().segment(2, n - 2));
    } 
   
    Vector cellCrossSection = (mixture.collision_data().totalCrossSection().segment(0, n) + mixture.collision_data().totalCrossSection().segment(1, n)) / 2.;

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
    const double Tg = m_workingConditions->gasTemperature();
    // reset by an assignment of a default-constructed Power object
    power = Power();
    elasticOperator.evaluatePower(grid(),eedf,Tg,power.elasticNet,power.elasticGain,power.elasticLoss);
    if (carOperator.get())
    {
        carOperator->evaluatePower(grid(),eedf,Tg,power.carNet,power.carGain,power.carLoss);
    }
    if (ionizationOperator.includeNonConservativeIonization || attachmentOperator.includeNonConservativeAttachment)
    {
        if (growthModelType == GrowthModelType::temporal)
        {
            // the calculation of field is identical to that in fieldOperator, except
            // for the usage of g_fieldTemporalGrowth instead of fieldOperator.g;
            // the former is based on a totalCS that has an additional term
            // growthFactor / std::sqrt(grid().getNode(i)).
            double field = 0., growthModel = 0.;
            for (Grid::Index k = 0; k < grid().nCells(); ++k)
            {
                field += eedf[k] * (g_fieldTemporalGrowth[k + 1] - g_fieldTemporalGrowth[k]);
                growthModel += eedf[k] * grid().getCell(k) * std::sqrt(grid().getCell(k));
            }
            power.field = SI::gamma * field;
            power.eDensGrowth = -CIEff * grid().du() * growthModel;
        }
        else if (growthModelType == GrowthModelType::spatial)
        {
            // first term 'field': same as in the case that there is no growth
            double field = 0.;
            fieldOperator.evaluatePower(grid(),eedf,field);
            // now calculate the additional terms
            double correction = 0., powerDiffusion = 0., powerMobility = 0.;
            Vector cellCrossSection(grid().nCells());
            for (Grid::Index k = 0; k < grid().nCells(); ++k)
            {
                correction -= eedf[k] * (g_fieldSpatialGrowth[k + 1] + g_fieldSpatialGrowth[k])/2;

                // Diffusion and Mobility contributions
                cellCrossSection[k] = .5 * (mixture.collision_data().totalCrossSection()[k] + mixture.collision_data().totalCrossSection()[k + 1]);
                powerDiffusion += grid().getCell(k) * grid().getCell(k) * eedf[k] / cellCrossSection[k];

                if (k > 0 && k < grid().nCells() - 1)
                {
                    powerMobility +=
                        grid().getCell(k) * grid().getCell(k) * (eedf[k + 1] - eedf[k - 1]) / cellCrossSection[k];
                }
            }

            /** \todo field and correction are both of the form 'eedf*g'.
             *  Check that the following is correct. That requires that g_E and g_fieldSpatialGrowth
             *  have different dimensions (factor energy).
             */
            // Note that field is calculated by fieldOperator.evaluatePower, which already
            // adds the factor SI::gamma.
            power.field = field + SI::gamma * grid().du() * correction;
            power.eDensGrowth = alphaRedEff * alphaRedEff * SI::gamma * grid().du() / 3. * powerDiffusion +
                                SI::gamma * alphaRedEff * (m_workingConditions->reducedFieldSI() / 6.) *
                                    (grid().getCell(0) * grid().getCell(0) * eedf[1] / cellCrossSection[0] -
                                    grid().getCell(grid().nCells() - 1) * grid().getCell(grid().nCells() - 1) *
                                        eedf[grid().nCells() - 2] / cellCrossSection[grid().nCells() - 1] +
                                    powerMobility);
        }
    }
    else
    {
        fieldOperator.evaluatePower(grid(),eedf,power.field);
    }

    if (eeOperator)
    {
        eeOperator->evaluatePower(grid(),eedf,power.electronElectron);
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
    double powerValues[13]{
        power.field,
        power.elasticGain, power.elasticLoss,
        power.carGain, power.carLoss,
        power.excitation.forward, power.excitation.backward,
        power.vibrational.forward, power.vibrational.backward,
        power.rotational.forward, power.rotational.backward,
        power.eDensGrowth,
        power.electronElectron
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
        + power.electronElectron;
    power.relativeBalance = std::abs(power.balance) / totalGain;
    power.reference = totalGain;
}

void ElectronKineticsBoltzmann::evaluateSwarmParameters()
{
    const Grid::Index n = grid().nCells();

    const bool nonConservative = (ionizationOperator.includeNonConservativeIonization || attachmentOperator.includeNonConservativeAttachment);

    Vector tCS(mixture.collision_data().totalCrossSection());

    if (growthModelType == GrowthModelType::temporal && nonConservative)
    {
        tCS.tail(grid().nCells()).array() += (CIEff/SI::gamma) / grid().getNodes().tail(n).cwiseSqrt().array();
    }
    if (grid().isUniform())
    {
        swarmParameters.redDiffCoeff = 2. / 3. * SI::gamma * grid().du() *
                                    grid().getCells().cwiseProduct(eedf).cwiseQuotient(tCS.head(n) + tCS.tail(n)).sum();
    } else
    {
        swarmParameters.redDiffCoeff = 2. / 3. * SI::gamma * grid().duCells().cwiseProduct(
                                    grid().getCells().cwiseProduct(eedf)).cwiseQuotient(tCS.head(n) + tCS.tail(n)).sum();
    }
    swarmParameters.redMobCoeff = -SI::gamma / 3. *
                                  grid().getNodes()
                                      .segment(1, n - 1)
                                      .cwiseProduct(eedf.tail(n - 1) - eedf.head(n - 1))
                                      .cwiseQuotient(tCS.segment(1, n - 1))
                                      .sum();

    if (growthModelType == GrowthModelType::spatial && nonConservative)
    {
        swarmParameters.driftVelocity = -swarmParameters.redDiffCoeff * alphaRedEff +
                                        swarmParameters.redMobCoeff * m_workingConditions->reducedFieldSI();
    }
    else
    {
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

    if (grid().isUniform())
    {
        swarmParameters.meanEnergy = grid().du() * (grid().getCells().array().pow(1.5) * eedf.array()).sum();
    } else 
    {
        swarmParameters.meanEnergy = (grid().duCells().array()*grid().getCells().array().pow(1.5)*eedf.array()).sum();
    }
    swarmParameters.characEnergy = swarmParameters.redDiffCoeff / swarmParameters.redMobCoeff;

    swarmParameters.Te = 2. / 3. * swarmParameters.meanEnergy;

    // TODO: is this correct? (simulations after the first will have a different value for Te).
    m_workingConditions->updateElectronTemperature(swarmParameters.Te);
}





ElectronKineticsPrescribed::ElectronKineticsPrescribed(const std::filesystem::path &basePath,const ElectronKineticsSetup &setup, WorkingConditions *workingConditions)
: ElectronKinetics(basePath, setup,workingConditions),
  shapeParameter(setup.shapeParameter)
{
    if (setup.ionizationOperatorType != IonizationOperatorType::conservative)
    {
        throw std::runtime_error("ionizationOperatorType must be 'conservative' for EEDF type 'Prescribed'.");
    }
    initialize();
}

ElectronKineticsPrescribed::ElectronKineticsPrescribed(const std::filesystem::path &basePath,const json_type &cnf, WorkingConditions *workingConditions)
: ElectronKinetics(basePath,cnf,workingConditions),
  shapeParameter(cnf.at("shapeParameter").get<double>())
{
    if (getIonizationOperatorType(cnf.at("ionizationOperatorType")) != IonizationOperatorType::conservative)
    {
        throw std::runtime_error("ionizationOperatorType must be 'conservative' for EEDF type 'Prescribed'.");
    }
    initialize();
}

void ElectronKineticsPrescribed::initialize()
{
    grid().updatedMaxEnergy.addListener(&ElectronKineticsPrescribed::evaluateMatrix, this);
    m_workingConditions->updatedReducedField.addListener(&ElectronKineticsPrescribed::evaluateFieldOperator, this);
    /// \todo Do this here? Not all parameters may have been set at this point.
    evaluateMatrix();
}

void ElectronKineticsPrescribed::evaluateFieldOperator()
{
    const double EoN = m_workingConditions->reducedFieldSI();
    const double WoN = m_workingConditions->reducedExcFreqSI();
    fieldOperator.evaluate(grid(),mixture.collision_data().totalCrossSection(),EoN,WoN);
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

    double totalGain = 0;
    double totalLoss = 0;

    /** \todo get rid of this; do 'for (double value : { ... this list ... })'
     *  in the line below.
     */
    double powerValues[13]{
        power.field,
        power.elasticGain, power.elasticLoss,
        power.carGain, power.carLoss,
        power.excitation.forward, power.excitation.backward,
        power.vibrational.forward, power.vibrational.backward,
        power.rotational.forward, power.rotational.backward,
        power.eDensGrowth,
        power.electronElectron
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
        + power.electronElectron;
    power.relativeBalance = std::abs(power.balance) / totalGain;
    power.reference = totalGain;
}

void ElectronKineticsPrescribed::evaluateSwarmParameters()
{
    const Grid::Index n = grid().nCells();

    Vector tCS(mixture.collision_data().totalCrossSection());

    swarmParameters.redDiffCoeff = 2. / 3. * SI::gamma * grid().du() *
                                   grid().getCells().cwiseProduct(eedf).cwiseQuotient(tCS.head(n) + tCS.tail(n)).sum();

    swarmParameters.redMobCoeff = -SI::gamma / 3. *
                                  grid().getNodes()
                                      .segment(1, n - 1)
                                      .cwiseProduct(eedf.tail(n - 1) - eedf.head(n - 1))
                                      .cwiseQuotient(tCS.segment(1, n - 1))
                                      .sum();

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

    swarmParameters.meanEnergy = grid().du() * (grid().getCells().array().pow(1.5) * eedf.array()).sum();

    swarmParameters.characEnergy = swarmParameters.redDiffCoeff / swarmParameters.redMobCoeff;

    swarmParameters.Te = 2. / 3. * swarmParameters.meanEnergy;

    // TODO: is this correct? (simulations after the first will have a different value for Te).
    m_workingConditions->updateElectronTemperature(swarmParameters.Te);
}

} // namespace loki
