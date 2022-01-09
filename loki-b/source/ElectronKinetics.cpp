//
// Created by daan on 13-5-19.
//

#include "LoKI-B/ElectronKinetics.h"
#include "LoKI-B/Constant.h"
#include <chrono>
#include <cmath>
#include <iomanip>

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

ElectronKinetics::ElectronKinetics(const ElectronKineticsSetup &setup, WorkingConditions *workingConditions)
    : workingConditions(workingConditions),
      grid(setup.numerics.energyGrid),
      mixture(&grid, setup, workingConditions),
      elasticMatrix(grid.nCells(), grid.nCells()),
      g_c(grid.getNodes().size()),
      attachmentConservativeMatrix(grid.nCells(), grid.nCells()),
      attachmentMatrix(grid.nCells(), grid.nCells()),
      fieldMatrix(grid.nCells(), grid.nCells()),
      g_E(grid.getNodes().size()),
      fieldMatrixSpatGrowth(grid.nCells(), grid.nCells()),
      ionSpatialGrowthD(grid.nCells(), grid.nCells()),
      ionSpatialGrowthU(grid.nCells(), grid.nCells()),
      fieldMatrixTempGrowth(grid.nCells(), grid.nCells()),
      ionTemporalGrowth(grid.nCells(), grid.nCells()),
      boltzmannMatrix(grid.nCells(), grid.nCells()),
      eedf(grid.nCells())
{
    this->eedfType = setup.eedfType;
    this->shapeParameter = setup.shapeParameter;
    this->mixingParameter = setup.numerics.nonLinearRoutines.mixingParameter;
    this->maxEedfRelError = setup.numerics.nonLinearRoutines.maxEedfRelError;
    this->maxPowerBalanceRelError = setup.numerics.maxPowerBalanceRelError;
    this->ionizationOperatorType = setup.ionizationOperatorType;
    this->growthModelType = setup.growthModelType;
    this->includeEECollisions = setup.includeEECollisions;

    initialize();
}

ElectronKinetics::ElectronKinetics(const json_type &cnf, WorkingConditions *workingConditions)
    : workingConditions(workingConditions),
    grid(cnf.at("numerics").at("energyGrid")),
    mixture(&grid, cnf, workingConditions),
    elasticMatrix(grid.nCells(), grid.nCells()),
    g_c(grid.nCells()),
    attachmentConservativeMatrix(grid.nCells(), grid.nCells()),
    attachmentMatrix(grid.nCells(), grid.nCells()),
    fieldMatrix(grid.nCells(), grid.nCells()),
    g_E(grid.getNodes().size()),
    fieldMatrixSpatGrowth(grid.nCells(), grid.nCells()),
    ionSpatialGrowthD(grid.nCells(), grid.nCells()),
    ionSpatialGrowthU(grid.nCells(), grid.nCells()),
    fieldMatrixTempGrowth(grid.nCells(), grid.nCells()),
    ionTemporalGrowth(grid.nCells(), grid.nCells()),
    boltzmannMatrix(grid.nCells(), grid.nCells()),
    eedf(grid.nCells())
{
    this->eedfType = getEedfType(cnf.at("eedfType"));
    this->shapeParameter = cnf.contains("shapeParameter") ? cnf.at("shapeParameter").get<unsigned>() : 0;
    this->mixingParameter = cnf.at("numerics").at("nonLinearRoutines").at("mixingParameter");
    this->maxEedfRelError = cnf.at("numerics").at("nonLinearRoutines").at("maxEedfRelError");
    this->maxPowerBalanceRelError = cnf.at("numerics").at("maxPowerBalanceRelError");
    this->ionizationOperatorType = getIonizationOperatorType(cnf.at("ionizationOperatorType"));
    this->growthModelType = getGrowthModelType(cnf.at("growthModelType"));
    this->includeEECollisions = cnf.at("includeEECollisions");

    initialize();
}

void ElectronKinetics::initialize()
{
    grid.updatedMaxEnergy.addListener(&ElectronKinetics::evaluateMatrix, this);
    workingConditions->updatedReducedField.addListener(&ElectronKinetics::evaluateFieldOperator, this);

    // this->plot("Total Elastic Cross Section N2", "Energy (eV)", "Cross Section (m^2)",
    // mixture.grid->getNodes(), mixture.collision_data().totalCrossSection());

    inelasticMatrix.setZero(grid.nCells(), grid.nCells());

    ionConservativeMatrix.setZero(grid.nCells(), grid.nCells());

    attachmentConservativeMatrix.setZero(grid.nCells(), grid.nCells());

    if (ionizationOperatorType != IonizationOperatorType::conservative &&
        mixture.collision_data().hasCollisions(CollisionType::ionization))
    {
        ionizationMatrix.setZero(grid.nCells(), grid.nCells());
    }

    A.setZero(grid.nCells());
    B.setZero(grid.nCells());

    boltzmannMatrix.setZero();

    // SPARSE INITIALIZATION
    std::vector<Eigen::Triplet<double>> tridiagPattern;
    tridiagPattern.reserve(3 * grid.nCells() - 2);

    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        if (k > 0)
            tridiagPattern.emplace_back(k - 1, k, 0.);

        tridiagPattern.emplace_back(k, k, 0.);

        if (k < grid.nCells() - 1)
            tridiagPattern.emplace_back(k + 1, k, 0.);
    }

    elasticMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());
    fieldMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());

    if (!mixture.CARGases().empty())
        CARMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());

    /** \todo Could the matrices that are guaranteed to be diagonal just be Eigen::DiagonalMatrix?
     *        That saves a lot of time initializing and makes accessing easier (no 'coeffRef').
     *        Then the pattern-code below can be removed, only a resize needed. Question:
     *        does Eigen optimize matrix addition and multiplication for such matrices? That
     *        is not immediately clear to me.
     *        Same for the other constructor, of course.
     */
    std::vector<Eigen::Triplet<double>> diagPattern;
    diagPattern.reserve(grid.nCells());

    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        diagPattern.emplace_back(k, k, 0.);
    }

    if (mixture.collision_data().hasCollisions(CollisionType::attachment))
    {
        attachmentMatrix.setFromTriplets(diagPattern.begin(), diagPattern.end());
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

    this->evaluateMatrix();
}

void ElectronKinetics::solveSingle()
{
    if (includeNonConservativeIonization || includeNonConservativeAttachment || includeEECollisions)
    {
        this->mixingDirectSolutions();
    }
    else
    {
        this->invertLinearMatrix();
    }
}

void ElectronKinetics::solveSmartGrid()
{
    assert (grid.smartGrid());
    const Grid::SmartGridParameters& smartGrid = *grid.smartGrid();
    solveSingle();
    double decades = calcDecades(eedf[0],eedf[grid.nCells()-1]);
    //std::cout << "decades: " << decades << ", uMax: " << grid.uMax() << std::endl;

    while (decades < smartGrid.minEedfDecay)
    {
        grid.updateMaxEnergy(grid.uMax() * (1 + smartGrid.updateFactor));
        solveSingle();
        decades = calcDecades(eedf[0],eedf[grid.nCells()-1]);
        //std::cout << "decades: " << decades << ", uMax: " << grid.uMax() << std::endl;
    }

    while (decades > smartGrid.maxEedfDecay)
    {
        grid.updateMaxEnergy(grid.uMax() / (1 + smartGrid.updateFactor));
        solveSingle();
        decades = calcDecades(eedf[0],eedf[grid.nCells()-1]);
        //std::cout << "decades: " << decades << ", uMax: " << grid.uMax() << std::endl;
    }
}

void ElectronKinetics::solveSmartGrid2()
{
    assert (grid.smartGrid());
    const Grid::SmartGridParameters& smartGrid = *grid.smartGrid();

    // 0. Set uM and uP equal to the present (initial) uMax, solve
    //    and calculate decades.
    double uM=grid.uMax();
    double uP=grid.uMax();
    solveSingle();
    double decades = calcDecades(eedf[0],eedf[grid.nCells()-1]);
    std::cout << "uMax = " << grid.uMax() << ", decades = " << decades << std::endl;
    // 1. Ensure that [decades(uM),decades(uP)] encloses [dM,dP],
    //    entirely or partially.
    if (decades<smartGrid.minEedfDecay)
    {
        // decades(uP) < dM: until decades(uP) >= dM, keep doubling uP,
        // set uMax=uP, solve and recalculate decades(uP)
        while (decades<smartGrid.minEedfDecay)
        {
            uP *= 2;
            grid.updateMaxEnergy(uP);
            solveSingle();
            decades = calcDecades(eedf[0],eedf[grid.nCells()-1]);
            std::cout << "uMax = " << grid.uMax() << ", decades = " << decades << std::endl;
        }
    }
    else if (decades>smartGrid.maxEedfDecay)
    {
        // decades(uM) > dP: until decades(uM) <= dP, keep dividing uM by 2,
        // set uMax=uM, solve and recalculate decades(uM)
        while (decades>smartGrid.maxEedfDecay)
        {
            uM /= 2;
            grid.updateMaxEnergy(uM);
            solveSingle();
            decades = calcDecades(eedf[0],eedf[grid.nCells()-1]);
            std::cout << "uMax = " << grid.uMax() << ", decades = " << decades << std::endl;
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
        grid.updateMaxEnergy((uP+uM)/2);
        solveSingle();
        decades = calcDecades(eedf[0],eedf[grid.nCells()-1]);
        std::cout << "uMax = " << grid.uMax() << ", decades = " << decades << std::endl;
        if (decades<smartGrid.minEedfDecay)
        {
            uM = grid.uMax();
        }
        else
        {
            uP = grid.uMax();
        }
    }
    std::cout << "Final uMax = " << grid.uMax() << ", decades = " << decades << std::endl;
}

void ElectronKinetics::solve()
{
    if (grid.smartGrid())
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
    mixture.collision_data().evaluateRateCoefficients(grid,eedf);
    evaluateSwarmParameters();
    evaluateFirstAnisotropy();

    obtainedNewEedf.emit(grid, eedf, *workingConditions, power, mixture.collision_data(), swarmParameters,
                         firstAnisotropy);

#ifdef LOKIB_CREATE_SPARSITY_PICTURE
    const std::string xpm_fname{"system_matrix.xpm"};
    std::cout << "Creating '" << xpm_fname << "'." << std::endl;
    writeXPM(boltzmannMatrix, xpm_fname);
#endif
}

const Grid *ElectronKinetics::getGrid()
{
    return &grid;
}

void ElectronKinetics::invertLinearMatrix()
{
    if (!mixture.CARGases().empty())
    {
        /// \todo Document all the scalings ('1e20') in this file. Are these really needed?
        boltzmannMatrix = 1.e20 * (elasticMatrix + fieldMatrix + CARMatrix + inelasticMatrix + ionConservativeMatrix +
                                   attachmentConservativeMatrix);
    }
    else
    {
        boltzmannMatrix = 1.e20 * (elasticMatrix + fieldMatrix + inelasticMatrix + ionConservativeMatrix +
                                   attachmentConservativeMatrix);
    }

    invertMatrix(boltzmannMatrix);
}

void ElectronKinetics::invertMatrix(Matrix &matrix)
{
//#define LOKIB_TIM_INVERT_MATRIX
#ifdef LOKIB_TIM_INVERT_MATRIX
    auto begin = std::chrono::high_resolution_clock::now();
#endif

    if (!hasSuperelastics)
    {
        eedf.setZero();
        eedf[0] = 1.;

        matrix.row(0) = grid.getCells().cwiseSqrt() * grid.du();
        /** We can also just make eedf[0]=1 by implementing the above and:
         *  matrix.row(0).setZero(); matrix(0,0)=1.0;
         *  The advantage is that the first row only has the diagonal,
         *  over-all the matrix has a better sparsity pattern.
         *  Then we can do the normalization afterwards.
         */

        LinAlg::hessenberg(matrix.data(), eedf.data(), grid.nCells());
    }
    else
    {
        /** \todo Explain the idea that is outlined below. Is it correct that p is not used
         *        after the call to LinAlg::hessenbergReductionPartialPiv? Is this icomplete?
         */
        // TODO: Find a way to distinguish when to use LU and when to use Hessenberg reduction.

        // HESSENBERG WITH PARTIAL PIVOTING
        //            auto *p = new uint32_t[grid.nCells()];
        //
        //            for (Grid::Index i = 0; i < grid.nCells(); ++i)
        //                p[i] = i;
        //
        //            LinAlg::hessenbergReductionPartialPiv(matrix.data(), &superElasticThresholds[0], p,
        //              grid.nCells(), superElasticThresholds.size());
        //
        //            eedf.setZero();
        //            eedf[0] = 1.;
        //
        //            matrix.row(0) = grid.getCells().cwiseSqrt() * grid.du();
        //
        //            LinAlg::hessenberg(matrix.data(), eedf.data(), grid.nCells());
        //
        //            delete[] p;

        // LU DECOMPOSITION
        Vector b = Vector::Zero(grid.nCells());
        b[0] = 1;

        matrix.row(0) = grid.getCells().cwiseSqrt() * grid.du();
        eedf = matrix.partialPivLu().solve(b);
    }

    /** \todo It seems that the normaization is superfluous, since the normalization condition
     *        is already part of the system (first row of A, first element of b). One could
     *        decide to change the first equation into eedf[0] = 1 and do the normalization
     *        afterwards. That prevents a fully populated first row of the system matrix
     *        (better sparsity pattern).
     */
    // std::cout << "NORM: " <<  eedf.dot(grid.getCells().cwiseSqrt() * grid.du()) << std::endl;
    eedf /= eedf.dot(grid.getCells().cwiseSqrt() * grid.du());

#ifdef LOKIB_TIM_INVERT_MATRIX
    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "Inverted matrix elapsed time = "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;
#endif
}

void ElectronKinetics::evaluateMatrix()
{
    mixture.collision_data().evaluateTotalAndElasticCS(grid);

    evaluateElasticOperator();

    evaluateFieldOperator();

    if (!mixture.CARGases().empty())
        evaluateCAROperator();

    evaluateInelasticOperators();

    if (mixture.collision_data().hasCollisions(CollisionType::ionization))
        evaluateIonizationOperator();

    if (mixture.collision_data().hasCollisions(CollisionType::attachment))
        evaluateAttachmentOperator();

    /** \todo see the comments about superElasticThresholds in the header file.
    // Sort and erase duplicates.
    std::sort(superElasticThresholds.begin(), superElasticThresholds.end());
    superElasticThresholds.erase(unique(superElasticThresholds.begin(), superElasticThresholds.end()),
                                 superElasticThresholds.end());
    */
}

void ElectronKinetics::evaluateElasticOperator()
{
    g_c = grid.getNodes().cwiseAbs2().cwiseProduct(mixture.collision_data().elasticCrossSection()) * 2;
    g_c[0] = 0.;
    g_c[g_c.size() - 1] = 0.;

    const double Tg = workingConditions->gasTemperature();
    const double c_el = Constant::kBeV * Tg;

    const double factor1 = (c_el / grid.du() + 0.5) / grid.du();
    const double factor2 = (c_el / grid.du() - 0.5) / grid.du();
    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        elasticMatrix.coeffRef(k, k) = -(g_c[k] * factor1 + g_c[k + 1] * factor2);

        if (k > 0)
            elasticMatrix.coeffRef(k, k - 1) = g_c[k] * factor2;

        if (k < grid.nCells() - 1)
            elasticMatrix.coeffRef(k, k + 1) = g_c[k + 1] * factor1;
    }
}

void ElectronKinetics::evaluateFieldOperator()
{
    const double EoN = workingConditions->reducedFieldSI();
    const double WoN = workingConditions->reducedExcFreqSI();
    const double me = Constant::electronMass;
    const double e = Constant::electronCharge;

    const Vector &cs = mixture.collision_data().totalCrossSection();

    // g_E gets its size in the constructor.
    g_E[0] = 0.;
    for (Grid::Index i=1; i!= g_E.size()-1; ++i)
    {
        g_E[i] = (EoN * EoN / 3) * grid.getNode(i) /
          (cs[i] + (me * WoN * WoN / (2 * e)) / (grid.getNode(i)*cs[i]));
    }
    g_E[g_E.size() - 1] = 0.;

    const double sqStep = grid.du() * grid.du();

    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        fieldMatrix.coeffRef(k, k) = -(g_E[k] + g_E[k + 1]) / sqStep;

        if (k > 0)
            fieldMatrix.coeffRef(k, k - 1) = g_E[k] / sqStep;

        if (k < grid.nCells() - 1)
            fieldMatrix.coeffRef(k, k + 1) = g_E[k + 1] / sqStep;
    }
}

/** \todo In the code all the G's are divided by N*sqrt(2*e/m_e),
 *  compared to the LoKI-B paper \cite Tejero2019, it seems.
 *  That explains why, in the code below, you see gas->fraction,
 *  whereas in the paper you see N_k. It would be good to have a
 *  document where the equations are written *exactly* as in the code.
 *  Also g is defined without the minus sign that appears in the definition
 *  in the paper. All in all, CARmatrix seems to be defined such that
 *  [CARmatrix]*[f] is an approximation of -(1/(N*sqrt(2*e/m_e))dG_CAR/du.
 */
void ElectronKinetics::evaluateCAROperator()
{
    /* When comparing this with Tejero2019, realize that in that paper,
     * equation 6c, the following symbols are used for a gas k:
     *
     * B_k: rotationalConstant
     * Q_{k,au}: quadruple moment in (atomic) units e*a_0^2. Here e is
     *   the elementary charge and a_0 the Bohr radius. (NOTE that the
     *   variable electricQuadrupoleMoment in the code is in SI units Cm^2.)
     * sigma_{0,k} = (8./15)*pi*Q_{k,au}^2*a_0^2, see \cite Tejero below
     * equation 6d, \cite Ridenti below equation 8b or Gerjuoy and Stein,
     * equation 20.
     *
     * For mixtures, the terms B_k*sigma_k in the expression for g_CAR must be
     * weighted with the molar fractions. The code first calculates this weighted
     * sum sigma0B, which gives
     *
     *   sum_k chi_k*B_k*sigma_k = (8./15)*pi*a_0^2* sum_k chi_k*B_k*Q_{k,au}^2
     *
     * NOTE: in the first public release of LoKI-B Q_{k,au} was used instead
     *       of Q_{k,au}^2. That will be fixed in version 2.0.0.
     */
    const double a02 = Constant::bohrRadius*Constant::bohrRadius;
    double sigma0B = 0.;
    for (const auto &gas : mixture.CARGases())
    {
        const double Qau = gas->electricQuadrupoleMoment/(Constant::electronCharge*a02);
        sigma0B += gas->fraction * Qau * Qau * gas->rotationalConstant;
    }
    sigma0B *= (8.*Constant::pi*a02/15.)*sigma0B;

    g_CAR = grid.getNodes() * (4. * sigma0B);
    // Boundary conditions. See Tejero2019 below equation 16b.
    g_CAR[0] = 0.;
    g_CAR[grid.nCells()] = 0.;

    const double Tg = workingConditions->gasTemperature();
    const double c_CAR = Constant::kBeV * Tg;

    const double factor1 = (c_CAR / grid.du() + 0.5) / grid.du();
    const double factor2 = (c_CAR / grid.du() - 0.5) / grid.du();
    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        CARMatrix.coeffRef(k, k) = -(g_CAR[k] * factor1 + g_CAR[k + 1] * factor2);

        if (k > 0)
            CARMatrix.coeffRef(k, k - 1) = g_CAR[k] * factor2;

        if (k < grid.nCells() - 1)
            CARMatrix.coeffRef(k, k + 1) = g_CAR[k + 1] * factor1;
    }
}

void ElectronKinetics::evaluateInelasticOperators()
{
    const Grid::Index cellNumber = grid.nCells();

    inelasticMatrix.setZero();

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (auto vecIndex : { CollisionType::excitation, CollisionType::vibrational, CollisionType::rotational})
        {
            for (const auto &collision : cd.collisions(vecIndex) )
            {
                const double threshold = collision->crossSection->threshold();

                if (threshold < grid.du() || threshold > grid.getNodes()[grid.nCells()])
                    continue;

                const double targetDensity = collision->getTarget()->delta();

                if (targetDensity != 0)
                {
                    const auto numThreshold = static_cast<Grid::Index>(std::floor(threshold / grid.du()));

                    Vector cellCrossSection(cellNumber);

                    for (Grid::Index i = 0; i < cellNumber; ++i)
                        cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

                    for (Grid::Index k = 0; k < cellNumber; ++k)
                    {
                        if (k < cellNumber - numThreshold)
                            inelasticMatrix(k, k + numThreshold) +=
                                targetDensity * grid.getCells()[k + numThreshold] * cellCrossSection[k + numThreshold];

                        inelasticMatrix(k, k) -= targetDensity * grid.getCells()[k] * cellCrossSection[k];
                    }

                    if (collision->isReverse())
                    {
                        const double swRatio = collision->getTarget()->statisticalWeight /
                                               collision->m_rhsHeavyStates[0]->statisticalWeight;
                        const double productDensity = collision->m_rhsHeavyStates[0]->delta();

                        if (productDensity == 0)
                            continue;

                        /** \todo see the comments about superElasticThresholds in the header file.
                        if (numThreshold > 1)
                            superElasticThresholds.emplace_back(numThreshold);
                        */

                        if (numThreshold != 1)
                            hasSuperelastics = true;

                        for (Grid::Index k = 0; k < cellNumber; ++k)
                        {
                            if (k >= numThreshold)
                                inelasticMatrix(k, k - numThreshold) +=
                                    swRatio * productDensity * grid.getCells()[k] * cellCrossSection[k];

                            if (k < cellNumber - numThreshold)
                                inelasticMatrix(k, k) -= swRatio * productDensity * grid.getCells()[k + numThreshold] *
                                                         cellCrossSection[k + numThreshold];
                        }
                    }
                }
            }
        }
    }
}

void ElectronKinetics::evaluateIonizationOperator()
{
    bool hasValidCollisions = false;

    ionConservativeMatrix.setZero();

    if (ionizationOperatorType != IonizationOperatorType::conservative)
        ionizationMatrix.setZero();

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (const auto &collision : cd.collisions(CollisionType::ionization))
        {
            const double threshold = collision->crossSection->threshold();

            if (threshold > grid.getNode(grid.nCells()))
                continue;

            hasValidCollisions = true;

            const double delta = collision->getTarget()->delta();
            const auto numThreshold = static_cast<Grid::Index>(std::floor(threshold / grid.du()));

            Vector cellCrossSection(grid.nCells());

            for (Grid::Index i = 0; i < grid.nCells(); ++i)
                cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

            switch (ionizationOperatorType)
            {
            case IonizationOperatorType::conservative:
                break;

            case IonizationOperatorType::oneTakesAll:
                for (Grid::Index k = 0; k < grid.nCells(); ++k)
                {
                    if (k < grid.nCells() - numThreshold)
                        ionizationMatrix(k, k + numThreshold) +=
                            delta * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                    const double term = delta * grid.getCell(k) * cellCrossSection(k);

                    ionizationMatrix(k, k) -= term;
                    ionizationMatrix(0, k) += term;
                }
                break;

            case IonizationOperatorType::equalSharing:
                for (Grid::Index k = 0; k < grid.nCells(); ++k)
                {
                    ionizationMatrix(k, k) -= delta * grid.getCell(k) * cellCrossSection[k];

                    if (k < (grid.nCells() - numThreshold) / 2)
                    {
                        const Grid::Index i = 2 * (k + 1) + numThreshold - 1;

                        ionizationMatrix(k, i) += 4 * delta * grid.getCell(i) * cellCrossSection(i);
                    }
                }
                break;
            case IonizationOperatorType::sdcs:
                double W = cd.OPBParameter();

                if (W < 0)
                    W = threshold;

                for (Grid::Index k = 0; k < grid.nCells(); ++k)
                {
                    const Grid::Index end = std::min(2 * (k + 1) + numThreshold, grid.nCells());

                    if (k > numThreshold)
                    {
                        const Grid::Index half = (k + 1 - numThreshold) / 2;
                        const double numerator = 1 / std::atan((grid.getCell(k) - threshold) / (2 * W));
                        double sum = 0.;

                        for (Grid::Index i = 0; i < half; ++i)
                            sum += numerator / (W + grid.getCell(i) * grid.getCell(i) / W);

                        ionizationMatrix(k, k) -= delta * grid.du() * grid.getCell(k) * cellCrossSection[k] * sum;
                    }

                    /** \todo If k + numThreshold + 1 >= grid.nCells(), the term is ignored.
                     *        Document (in the document, not necessarily here) what are the
                     *        consequences of that.
                     */
                    if (k + numThreshold + 1 < grid.nCells())
                    {
                        for (Grid::Index i = k + numThreshold + 1; i < end; ++i)
                        {
                            ionizationMatrix(k, i) += delta * grid.du() * grid.getCell(i) * cellCrossSection[i] /
                                                      (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                       (W + std::pow(grid.getCell(i - k - numThreshold - 1), 2) / W));
                        }
                    }

                    /** \todo The following comment needs to be sorted out (possible index errors).
                     *  \todo Document what is done here (algorithm) and how it is implemented.
                     */

                    /** \todo This last section might need some adjustments because of indexing
                     *  differences between Matlab and C++ (since indexes are multiplied here).
                     */

                    for (Grid::Index i = 2 * (k + 1) + numThreshold - 1; i < grid.nCells(); ++i)
                    {
                        ionizationMatrix(k, i) += delta * grid.du() * grid.getCell(i) * cellCrossSection[i] /
                                                  (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                   (W + std::pow(grid.getCell(k), 2) / W));
                    }
                }
                break;
            }

            // Evaluation of the conservative ionization operator

            if (numThreshold == 0)
                continue;

            for (Grid::Index k = 0; k < grid.nCells(); ++k)
            {
                if (k < grid.nCells() - numThreshold)
                    ionConservativeMatrix(k, k + numThreshold) +=
                        delta * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                ionConservativeMatrix(k, k) -= delta * grid.getCell(k) * cellCrossSection[k];
            }

            if (ionizationOperatorType != IonizationOperatorType::conservative && hasValidCollisions)
                includeNonConservativeIonization = true;
        }
    }
}

void ElectronKinetics::evaluateAttachmentOperator()
{
    attachmentMatrix.setZero();
    attachmentConservativeMatrix.setZero();

    const Grid::Index cellNumber = grid.nCells();

    for (const auto &cd : mixture.collision_data().data_per_gas())
    {
        for (const auto &collision : cd.collisions(CollisionType::attachment))
        {
            const double threshold = collision->crossSection->threshold();

            if (threshold > grid.getNode(cellNumber))
                continue;

            /* This should definitely not be in this (double) loop. Is this a constructor task?
             * Can this just be replaced with 'cd.collisions(CollisionType::attachment).size()'
             * in (other) places where this is now used?
             * Answer: no, this depends on uMax(), which may change for a smart grid. But it
             * could be done as an action when that changes.
             */
            /** \bug This should be reset to false at the beginning of this function
             *       because the results may change when uMax is changed.
             */
            includeNonConservativeAttachment = true;

            const auto numThreshold = static_cast<Grid::Index>(std::floor(threshold / grid.du()));

            /** \todo Eliminate the cellCrossSection vector? This can be calculated on the fly
             *        in the two places where it is needed (one if the merger below can be done).
             */
            Vector cellCrossSection(cellNumber);

            const double targetDensity = collision->getTarget()->delta();

            /// \todo Merge with the subsequent k-loop.
            for (Grid::Index i = 0; i < cellNumber; ++i)
                cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

            for (Grid::Index k = 0; k < cellNumber; ++k)
                attachmentMatrix.coeffRef(k, k) -= targetDensity * grid.getCell(k) * cellCrossSection[k];

            if (numThreshold == 0)
                continue;

            /** Can we also merge with this loop?
             *  It does not seem problematic to have an 'if (numThreshold)' in the loop,
             *  given the amount of work done.
             */
            for (Grid::Index k = 0; k < cellNumber; ++k)
            {
                if (k < cellNumber - numThreshold)
                    attachmentConservativeMatrix(k, k + numThreshold) +=
                        targetDensity * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                attachmentConservativeMatrix(k, k) -= targetDensity * grid.getCell(k) * cellCrossSection[k];
            }
        }
    }
}

void ElectronKinetics::mixingDirectSolutions()
{
    invertLinearMatrix();

    const Grid::Index numCells = grid.nCells();

    //        solveEEColl();
    //        return;

    // Declare function pointer
    void (ElectronKinetics::*growthFunc)() = nullptr;

    const bool includeGrowthModel = includeNonConservativeAttachment || includeNonConservativeIonization;

    if (includeGrowthModel)
    {
        switch (growthModelType)
        {
        case GrowthModelType::spatial:
            ionSpatialGrowthD.setZero();
            ionSpatialGrowthU.setZero();
            fieldMatrixSpatGrowth.setZero();

            growthFunc = &ElectronKinetics::solveSpatialGrowthMatrix;
            break;

        case GrowthModelType::temporal:
            ionTemporalGrowth.setZero();
            fieldMatrixTempGrowth.setZero();

            growthFunc = &ElectronKinetics::solveTemporalGrowthMatrix;
            break;
        }
    }

    if (includeEECollisions)
    {
        alphaEE = 0.;
        BAee.setZero(numCells, numCells);
        A.setZero();
        B.setZero();

        if (!includeGrowthModel)
        {
            Log<Message>::Notify("Starting e-e collision routine.");
            solveEEColl();
        }
        else
        {
            uint32_t globalIter = 0;
            const uint32_t maxGlobalIter = 20;

            Vector eedfOld;

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
    }
    else
    {
        (this->*growthFunc)();
    }
}

void ElectronKinetics::solveSpatialGrowthMatrix()
{
    const double EoN = workingConditions->reducedFieldSI();

    Vector cellTotalCrossSection(grid.nCells());

    for (Grid::Index i = 0; i < grid.nCells(); ++i)
        cellTotalCrossSection[i] = .5 * (mixture.collision_data().totalCrossSection()[i] + mixture.collision_data().totalCrossSection()[i + 1]);

    if (!mixture.CARGases().empty())
    {
        boltzmannMatrix = 1.e20 * (elasticMatrix + fieldMatrix + CARMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
    }
    else
    {
        boltzmannMatrix = 1.e20 * (elasticMatrix + fieldMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
    }

    Vector baseDiag(grid.nCells()), baseSubDiag(grid.nCells()), baseSupDiag(grid.nCells());

    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        baseDiag[k] = boltzmannMatrix(k, k);

        if (k > 0)
            baseSubDiag[k] = boltzmannMatrix(k, k - 1);

        if (k < grid.nCells() - 1)
            baseSupDiag[k] = boltzmannMatrix(k, k + 1);
    }

    if (includeEECollisions)
    {
        A = alphaEE / grid.du() * (BAee.transpose() * eedf);
        B = alphaEE / grid.du() * (BAee * eedf);
    }

    /** \todo The name is incorrect. This is not the integrand since it already includes
     *        du and does not yet have the factor f(u).
     */
    Vector integrandCI = (SI::gamma*grid.du()) * Vector::Ones(grid.nCells()).transpose() *
                         (ionizationMatrix + attachmentMatrix);

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

    // diffusion and mobility components of the spatial growth terms
    // can be removed since this is already done in the directMixing function
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
     *  if you have enough grid points).
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
    const Vector D0 = grid.getCells().array() / (3. * cellTotalCrossSection).array();
    /** \todo Document/explain which equation is discretized here.
     */
    Vector U0sup(grid.nCells());
    Vector U0inf(grid.nCells());
    U0sup[0] = 0.;
    U0inf[grid.nCells() - 1] = 0.;
    for (Grid::Index j = 0; j < grid.nCells(); ++j)
    {
        if (j != 0)
            U0sup[j] = EoN / (2. * grid.du()) * D0[j - 1];

        if (j != grid.nCells() - 1)
            U0inf[j] = -EoN / (2. * grid.du()) * D0[j + 1];
    }
    const Vector U0 = U0sup + U0inf;

    // This is 33a from \cite Manual_1_0_0
    double ND  =   SI::gamma * grid.du() * D0.dot(eedf);
    /* This is 33b from \cite Manual_1_0_0, multiplied with E/N.
     * Note that the factor E/N is part of U0sup, U0inf.
     */
    double muE = - SI::gamma * grid.du() * U0.dot(eedf);

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

    const Vector g_fieldSpatialBase = (EoN / 6) * grid.getNodes().array() / mixture.collision_data().totalCrossSection().array();

    while (!hasConverged)
    {
        g_fieldSpatialGrowth = alphaRedEffNew * g_fieldSpatialBase;
        g_fieldSpatialGrowth[0] = 0.;
        g_fieldSpatialGrowth[grid.nCells()] = 0.;

        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            fieldMatrixSpatGrowth.coeffRef(k, k) = (g_fieldSpatialGrowth[k + 1] - g_fieldSpatialGrowth[k]) / grid.du();

            if (k > 0)
                fieldMatrixSpatGrowth.coeffRef(k, k - 1) = -g_fieldSpatialGrowth[k] / grid.du();

            if (k < grid.nCells() - 1)
                fieldMatrixSpatGrowth.coeffRef(k, k + 1) = g_fieldSpatialGrowth[k + 1] / grid.du();
        }

        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            ionSpatialGrowthD.coeffRef(k, k) = alphaRedEffNew * alphaRedEffNew * D0[k];

            if (k > 0)
                ionSpatialGrowthU.coeffRef(k, k - 1) = alphaRedEffNew * U0inf[k - 1];

            if (k < grid.nCells() - 1)
                ionSpatialGrowthU.coeffRef(k, k + 1) = alphaRedEffNew * U0sup[k + 1];
        }

        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            boltzmannMatrix(k, k) = baseDiag[k] + 1.e20 * (fieldMatrixSpatGrowth.coeff(k, k) +
                                                           ionSpatialGrowthD.coeff(k, k) - (A[k] + B[k]));

            if (k > 0)
                boltzmannMatrix(k, k - 1) = baseSubDiag[k] + 1.e20 * (fieldMatrixSpatGrowth.coeff(k, k - 1) +
                                                                      ionSpatialGrowthU.coeff(k, k - 1) + A[k - 1]);

            if (k < grid.nCells() - 1)
                boltzmannMatrix(k, k + 1) = baseSupDiag[k] + 1.e20 * (fieldMatrixSpatGrowth.coeff(k, k + 1) +
                                                                      ionSpatialGrowthU.coeff(k, k + 1) + B[k + 1]);
        }

        Vector eedfNew = eedf;

        invertMatrix(boltzmannMatrix);

        CIEffOld = CIEffNew;
        CIEffNew = eedf.dot(integrandCI);

        CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

        ND  =   SI::gamma * grid.du() * D0.dot(eedf);
        muE = - SI::gamma * grid.du() * U0.dot(eedf);

        /** \todo See the notes just above this loop for a note about the discontinuity.
         */
        const double discriminant = muE * muE - 4 * CIEffNew * ND;

        alphaRedEffOld = alphaRedEffNew;
        alphaRedEffNew = (discriminant < 0) ? CIEffNew / muE : (muE - std::sqrt(discriminant)) / (2 * ND);

        alphaRedEffNew = mixingParameter * alphaRedEffNew + (1 - mixingParameter) * alphaRedEffOld;

        if (((alphaRedEffNew == 0 || std::abs(alphaRedEffNew - alphaRedEffOld) / alphaRedEffOld < 1.e-10) &&
             maxRelDiff(eedfNew,eedf) < maxEedfRelError) ||
            iter > 150)
        {
            hasConverged = true;

            /** There is no maximum number of iterations in case includeEECollisions==true.
             *  Is there a reason for that? This is not very safe.
             */
            if (iter > 150 && !includeEECollisions)
                Log<Message>::Warning("Iterative spatial growth scheme did not converge.");
        }
        ++iter;
    }

    std::cerr << "Spatial growth routine converged in: " << iter << " iterations.\n";

    alphaRedEff = alphaRedEffOld;
    CIEff = CIEffOld;
}

void ElectronKinetics::solveTemporalGrowthMatrix()
{
    const double e = Constant::electronCharge;
    const double m = Constant::electronMass;
    const double EoN = workingConditions->reducedFieldSI();
    const double WoN = workingConditions->reducedExcFreqSI();

    if (!mixture.CARGases().empty())
    {
        boltzmannMatrix = 1.e20 * (elasticMatrix + CARMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
    }
    else
    {
        boltzmannMatrix = 1.e20 * (elasticMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
    }

    Vector baseDiag(grid.nCells()), baseSubDiag(grid.nCells()), baseSupDiag(grid.nCells());

    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        baseDiag[k] = boltzmannMatrix(k, k);

        if (k > 0)
            baseSubDiag[k] = boltzmannMatrix(k, k - 1);

        if (k < grid.nCells() - 1)
            baseSupDiag[k] = boltzmannMatrix(k, k + 1);
    }

    if (includeEECollisions)
    {
        A = alphaEE / grid.du() * (BAee.transpose() * eedf);
        B = alphaEE / grid.du() * (BAee * eedf);
    }

    const Vector integrandCI = (SI::gamma * grid.du())
                    * Vector::Ones(grid.nCells()).transpose()
                    * (ionizationMatrix + attachmentMatrix);
    double CIEffNew = eedf.dot(integrandCI);
    double CIEffOld = CIEffNew / 3.;
    CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

    Vector eedfNew(grid.nCells());

    bool hasConverged = false;
    uint32_t iter = 0;

    while (!hasConverged)
    {
 //       Log<Message>::Notify("Iteration ", iter);

        const long double growthFactor = CIEffNew / SI::gamma;

        g_fieldTemporalGrowth.resize(grid.getNodes().size());
        g_fieldTemporalGrowth[0] = 0.;
        for (Grid::Index i=1; i!= g_fieldTemporalGrowth.size()-1; ++i)
        {
            const double totalCSI = mixture.collision_data().totalCrossSection()[i] + growthFactor / std::sqrt(grid.getNode(i));
            g_fieldTemporalGrowth[i] =
                (EoN * EoN / 3) * grid.getNode(i) /
                (totalCSI + (m * WoN * WoN / (2 * e)) / (grid.getNode(i)*totalCSI));
        }
        g_fieldTemporalGrowth[g_fieldTemporalGrowth.size() - 1] = 0.;

        const double sqrStep = grid.du() * grid.du();

        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            fieldMatrixTempGrowth.coeffRef(k, k) = -(g_fieldTemporalGrowth[k] + g_fieldTemporalGrowth[k + 1]) / sqrStep;

            if (k > 0)
                fieldMatrixTempGrowth.coeffRef(k, k - 1) = g_fieldTemporalGrowth[k] / sqrStep;

            if (k < grid.nCells() - 1)
                fieldMatrixTempGrowth.coeffRef(k, k + 1) = g_fieldTemporalGrowth[k + 1] / sqrStep;

            ionTemporalGrowth.coeffRef(k, k) = -growthFactor * std::sqrt(grid.getCell(k));
        }

        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            boltzmannMatrix(k, k) = baseDiag[k] + 1.e20 * (fieldMatrixTempGrowth.coeff(k, k) +
                                                           ionTemporalGrowth.coeff(k, k) - (A[k] + B[k]));

            if (k > 0)
                boltzmannMatrix(k, k - 1) = baseSubDiag[k] + 1.e20 * (fieldMatrixTempGrowth.coeff(k, k - 1) + A[k - 1]);

            if (k < grid.nCells() - 1)
                boltzmannMatrix(k, k + 1) = baseSupDiag[k] + 1.e20 * (fieldMatrixTempGrowth.coeff(k, k + 1) + B[k + 1]);
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

            if (iter > 150 && !includeEECollisions)
                Log<Message>::Warning("Iterative temporal growth scheme did not converge.");
        }

        ++iter;
    }

 //   std::cerr << "Temporal growth routine converged in: " << iter << " iterations.\n";

    CIEff = CIEffOld;
}

void ElectronKinetics::solveEEColl()
{
    const double e = Constant::electronCharge;
    const double e0 = Constant::vacuumPermittivity;
    const double ne = workingConditions->electronDensity();
    const double n0 = workingConditions->gasDensity();

    // Splitting all possible options for the best performance.

    if (includeNonConservativeIonization || includeNonConservativeAttachment)
    {
        if (growthModelType == GrowthModelType::spatial)
        {
            if (mixture.CARGases().empty())
            {
                boltzmannMatrix = 1.e20 * (ionizationMatrix + attachmentMatrix + elasticMatrix + inelasticMatrix +
                                           fieldMatrix + ionSpatialGrowthD + ionSpatialGrowthU + fieldMatrixSpatGrowth);
            }
            else
            {
                boltzmannMatrix =
                    1.e20 * (ionizationMatrix + attachmentMatrix + elasticMatrix + inelasticMatrix + CARMatrix +
                             fieldMatrix + ionSpatialGrowthD + ionSpatialGrowthU + fieldMatrixSpatGrowth);
            }
        }
        else if (growthModelType == GrowthModelType::temporal)
        {
            if (mixture.CARGases().empty())
            {
                boltzmannMatrix = 1.e20 * (ionizationMatrix + attachmentMatrix + elasticMatrix + inelasticMatrix +
                                           ionTemporalGrowth + fieldMatrixTempGrowth);
            }
            else
            {
                boltzmannMatrix = 1.e20 * (ionizationMatrix + attachmentMatrix + elasticMatrix + inelasticMatrix +
                                           CARMatrix + ionTemporalGrowth + fieldMatrixTempGrowth);
            }
        }
    }
    else
    {
        if (mixture.CARGases().empty())
        {
            boltzmannMatrix = 1.e20 * (ionConservativeMatrix + attachmentConservativeMatrix + elasticMatrix +
                                       inelasticMatrix + fieldMatrix);
        }
        else
        {
            boltzmannMatrix = 1.e20 * (ionConservativeMatrix + attachmentConservativeMatrix + elasticMatrix +
                                       inelasticMatrix + fieldMatrix + CARMatrix);
        }
    }

    // Storing the initial diagonals of the matrix in three separate vectors.
    // This allows us to skip the usage of 'matrixAux', saving a good amount of
    // memory for bigger matrices.

    Vector baseDiag(grid.nCells()), baseSubDiag(grid.nCells()), baseSupDiag(grid.nCells());

    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        baseDiag[k] = boltzmannMatrix(k, k);

        if (k > 0)
            baseSubDiag[k] = boltzmannMatrix(k, k - 1);

        if (k < grid.nCells() - 1)
            baseSupDiag[k] = boltzmannMatrix(k, k + 1);
    }

    BAee.setZero(grid.nCells(), grid.nCells());

    const Vector cellsThreeOverTwo = grid.getCells().cwiseProduct(grid.getCells().cwiseSqrt()),
                 energyArray = -(grid.du() / 2.) * grid.getCells().cwiseSqrt() + (2. / 3.) * cellsThreeOverTwo;

    for (Grid::Index j = 0; j < grid.nCells() - 1; ++j)
    {

        for (Grid::Index i = 1; i <= j; ++i)
            BAee(i, j) = energyArray[i];

        const double value = 2. / 3. * std::pow(grid.getNode(j + 1), 1.5);

        for (Grid::Index i = j + 1; i < grid.nCells(); ++i)
            BAee(i, j) = value;
    }

    // detailed balance condition

    for (Grid::Index j = 0; j < grid.nCells() - 1; ++j)
    {
        for (Grid::Index i = 1; i < grid.nCells(); ++i)
        {
            BAee(i, j) = std::sqrt(BAee(i, j) * BAee(j + 1, i - 1));
        }
    }

    /** \todo These four declarations can be moved to the beginning of the
     *        while loop. Then the calculation does not have to repeated
     *        at the end of that loop (and the decls. can be made constant).
     */
    double meanEnergy = grid.du() * cellsThreeOverTwo.dot(eedf);
    double Te = 2. / 3. * meanEnergy;
    double logC = std::log(12 * Constant::pi * std::pow(e0 * Te / e, 1.5) / std::sqrt(ne));
    double alpha = (ne / n0) * (e * e / (8 * Constant::pi * e0 * e0)) * logC;

    double ratioNew = 0.;
    Vector eedfNew = eedf;

    bool hasConverged = false;
    uint32_t iter = 0;

    /** Mee is ignored in solveEEColl. How does that result in differences with
     *  respect to the MATLAB version.
     */
    /// In this implementation we completely skip the Mee matrix, saving both memory and time.
    // Vector MeeDiag(grid.nCells()), MeeSub(grid.nCells()), MeeSup(grid.nCells());

    while (!hasConverged)
    {
        A = (alpha / grid.du()) * (BAee.transpose() * eedf);
        B = (alpha / grid.du()) * (BAee * eedf);

        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            boltzmannMatrix(k, k) = baseDiag[k] - 1.e20 * (A[k] + B[k]);

            if (k > 0)
                boltzmannMatrix(k, k - 1) = baseSubDiag[k] + 1.e20 * A[k - 1];

            if (k < grid.nCells() - 1)
                boltzmannMatrix(k, k + 1) = baseSupDiag[k] + 1.e20 * B[k + 1];
        }

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
        else if (iter == 300 && !(includeNonConservativeAttachment || includeNonConservativeIonization))
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
            for (Grid::Index i = 0; i < grid.nCells(); ++i)
            {
                if (eedf[i] < 0)
                {
                    eedf[i] = std::abs(eedf[i]);
                }
            }
        }

        meanEnergy = grid.du() * cellsThreeOverTwo.dot(eedf);
        Te = 2. / 3. * meanEnergy;
        logC = std::log(12 * Constant::pi * std::pow(e0 * Te / e, 1.5) / std::sqrt(ne));
        alpha = (ne / n0) * (e * e / (8 * Constant::pi * e0 * e0)) * logC;

        iter++;
    }

    alphaEE = alpha;

    std::cerr << "e-e routine converged in: " << iter << " iterations.\n";
}

void ElectronKinetics::evaluatePower()
{
    const double kTg = Constant::kBeV * workingConditions->gasTemperature();
    const double auxHigh = kTg + grid.du() * .5; // aux1
    const double auxLow  = kTg - grid.du() * .5; // aux2

    // reset by an assignment to a default-constructed Power object
    power = Power();

    double elasticNet = 0., elasticGain = 0.;
    for (Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        elasticNet += eedf[k] * (g_c[k + 1] * auxLow - g_c[k] * auxHigh);
        elasticGain += eedf[k] * (g_c[k + 1] - g_c[k]);
    }
    power.elasticNet = SI::gamma * elasticNet;
    power.elasticGain = SI::gamma * kTg * elasticGain;
    power.elasticLoss = power.elasticNet - power.elasticGain;

    if (!mixture.CARGases().empty())
    {
        double carNet = 0., carGain = 0.;
        for (Grid::Index k = 0; k < grid.nCells() - 1; ++k)
        {
            carNet += eedf[k] * (g_CAR[k + 1] * auxLow - g_CAR[k] * auxHigh);
            carGain += eedf[k] * (g_CAR[k + 1] - g_CAR[k]);
        }
        power.carNet = SI::gamma * carNet;
        power.carGain = SI::gamma * kTg * carGain;
        power.carLoss = power.carNet - power.carGain;
    }

    if (includeNonConservativeIonization || includeNonConservativeAttachment)
    {
        if (growthModelType == GrowthModelType::temporal)
        {
            double field = 0., growthModel = 0.;
            for (Grid::Index k = 0; k < grid.nCells(); ++k)
            {
                field += eedf[k] * (g_fieldTemporalGrowth[k + 1] - g_fieldTemporalGrowth[k]);
                growthModel += eedf[k] * grid.getCell(k) * std::sqrt(grid.getCell(k));
            }
            power.field = SI::gamma * field;
            power.eDensGrowth = -CIEff * grid.du() * growthModel;
        }
        else if (growthModelType == GrowthModelType::spatial)
        {
            double field = 0., correction = 0., powerDiffusion = 0., powerMobility = 0.;
            Vector cellCrossSection(grid.nCells());

            for (Grid::Index k = 0; k < grid.nCells(); ++k)
            {
                /// \todo check that it is correct that g_E (local) is used in a spatial growth simulation)
                field += eedf[k] * (g_E[k + 1] - g_E[k]);
                correction -= eedf[k] * (g_fieldSpatialGrowth[k + 1] + g_fieldSpatialGrowth[k]);

                // Diffusion and Mobility contributions
                cellCrossSection[k] = .5 * (mixture.collision_data().totalCrossSection()[k] + mixture.collision_data().totalCrossSection()[k + 1]);
                powerDiffusion += grid.getCell(k) * grid.getCell(k) * eedf[k] / cellCrossSection[k];

                if (k > 0 && k < grid.nCells() - 1)
                {
                    powerMobility +=
                        grid.getCell(k) * grid.getCell(k) * (eedf[k + 1] - eedf[k - 1]) / cellCrossSection[k];
                }
            }

            /** \todo field and correction are both of the form 'eedf*g'.
             *  Check that the following is correct. That requires that g_E and g_fieldSpatialGrowth
             *  have different dimensions (factor energy).
             */
            power.field = SI::gamma * (field + grid.du() * correction);
            power.eDensGrowth = alphaRedEff * alphaRedEff * SI::gamma * grid.du() / 3. * powerDiffusion +
                                SI::gamma * alphaRedEff * (workingConditions->reducedFieldSI() / 6.) *
                                    (grid.getCell(0) * grid.getCell(0) * eedf[1] / cellCrossSection[0] -
                                     grid.getCell(grid.nCells() - 1) * grid.getCell(grid.nCells() - 1) *
                                         eedf[grid.nCells() - 2] / cellCrossSection[grid.nCells() - 1] +
                                     powerMobility);
        }
    }
    else
    {
        double field = 0;
        for (Grid::Index k = 0; k < grid.nCells(); ++k)
        {
            field += eedf[k] * (g_E[k + 1] - g_E[k]);
        }
        power.field = SI::gamma * field;
    }

    if (includeEECollisions)
    {
        power.electronElectron = (-SI::gamma * grid.du() * grid.du()) * (A - B).dot(eedf);
    }

    // Evaluate power absorbed per electron at unit gas density due to in- and superelastic collisions.
    for (auto &cd : mixture.collision_data().data_per_gas())
    {
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

    double totalGain = 0., totalLoss = 0.;

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

void ElectronKinetics::evaluateSwarmParameters()
{
    const Grid::Index n = grid.nCells();

    const bool nonConservative = (includeNonConservativeIonization || includeNonConservativeAttachment);

    Vector tCS(mixture.collision_data().totalCrossSection());

    if (growthModelType == GrowthModelType::temporal && nonConservative)
    {
        tCS.tail(grid.nCells()).array() += (CIEff/SI::gamma) / grid.getNodes().tail(n).cwiseSqrt().array();
    }

    swarmParameters.redDiffCoeff = 2. / 3. * SI::gamma * grid.du() *
                                   grid.getCells().cwiseProduct(eedf).cwiseQuotient(tCS.head(n) + tCS.tail(n)).sum();

    swarmParameters.redMobCoeff = -SI::gamma / 3. *
                                  grid.getNodes()
                                      .segment(1, n - 1)
                                      .cwiseProduct(eedf.tail(n - 1) - eedf.head(n - 1))
                                      .cwiseQuotient(tCS.segment(1, n - 1))
                                      .sum();

    if (growthModelType == GrowthModelType::spatial && nonConservative)
    {
        swarmParameters.driftVelocity = -swarmParameters.redDiffCoeff * alphaRedEff +
                                        swarmParameters.redMobCoeff * workingConditions->reducedFieldSI();
    }
    else
    {
        swarmParameters.driftVelocity = swarmParameters.redMobCoeff * workingConditions->reducedFieldSI();
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

    swarmParameters.meanEnergy = grid.du() * (grid.getCells().array().pow(1.5) * eedf.array()).sum();

    swarmParameters.characEnergy = swarmParameters.redDiffCoeff / swarmParameters.redMobCoeff;

    swarmParameters.Te = 2. / 3. * swarmParameters.meanEnergy;

    // TODO: is this correct? (simulations after the first will have a different value for Te).
    workingConditions->updateElectronTemperature(swarmParameters.Te);
}

void ElectronKinetics::evaluateFirstAnisotropy()
{
    firstAnisotropy.setZero(grid.nCells());

    const double e = Constant::electronCharge;
    const double me = Constant::electronMass;
    const double EoN = workingConditions->reducedFieldSI();
    const double WoN = workingConditions->reducedExcFreqSI();
    const Grid::Index n = grid.nCells();

    // 1. First fill firstAnisotropy with df/du.
    firstAnisotropy[0] = (eedf[1] - eedf[0]) / grid.du();
    firstAnisotropy[n - 1] = (eedf[n - 1] - eedf[n - 2]) / grid.du();
    firstAnisotropy.segment(1, n - 2) = (eedf.segment(2, n - 2) - eedf.segment(0, n - 2)) / (2 * grid.du());

    Vector cellCrossSection = (mixture.collision_data().totalCrossSection().segment(0, n) + mixture.collision_data().totalCrossSection().segment(1, n)) / 2.;

    if (includeNonConservativeIonization || includeNonConservativeAttachment)
    {
        if (growthModelType == GrowthModelType::temporal)
        {
            // The DC case. This implements Tejero2019, equation 3b:
            // f1 = -zeta(E/N)(1/Omega_PT)(df/du).

            // Calculate Omega_c. This is the first term that appears in
            // Tejero2019 equation 5a. It is defined in the text below 5b.
            cellCrossSection =
                cellCrossSection.array() + CIEff / (SI::gamma*grid.getCells().cwiseSqrt()).array();
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
                // Omega_PT_i = Omega_c_i + (me/2))(omega/N)^2 / (e*u*Omega_c_i)
                // (Note: that e*u is the energy in SI units.)
                firstAnisotropy = -EoN * std::sqrt(2.) * firstAnisotropy.array() /
                                  (cellCrossSection.array() +
                                   WoN * WoN * me / (2. * e * grid.getCells().array() * cellCrossSection.array()));
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
        firstAnisotropy =
            -EoN * std::sqrt(2.) * firstAnisotropy.array() /
            (cellCrossSection.array() + WoN * WoN * me / (2. * e * grid.getCells().array() * cellCrossSection.array()));
    }
}
} // namespace loki

// Code to print a matrix
//            for (Grid::Index i = 0; i < grid.nCells(); ++i) {
//                for (Grid::Index j = 0; j < grid.nCells(); ++j) {
//                    printf("%.16e", matrix(i, j));
//
//                    if (j < grid.nCells() - 1) {
//                        printf("\t");
//                    }
//                }
//
//                printf("\n");
//            }
