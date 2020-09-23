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
    : workingConditions(workingConditions), grid(setup.numerics.energyGrid), mixture(&grid, setup, workingConditions),
      attachmentConservativeMatrix(grid.cellNumber, grid.cellNumber), boltzmannMatrix(grid.cellNumber, grid.cellNumber),
      elasticMatrix(grid.cellNumber, grid.cellNumber), fieldMatrix(grid.cellNumber, grid.cellNumber),
      attachmentMatrix(grid.cellNumber, grid.cellNumber), ionSpatialGrowthD(grid.cellNumber, grid.cellNumber),
      ionSpatialGrowthU(grid.cellNumber, grid.cellNumber), fieldMatrixSpatGrowth(grid.cellNumber, grid.cellNumber),
      fieldMatrixTempGrowth(grid.cellNumber, grid.cellNumber), ionTemporalGrowth(grid.cellNumber, grid.cellNumber),
      g_c(grid.cellNumber), eedf(grid.cellNumber)
{
    grid.updatedMaxEnergy2.addListener(&ElectronKinetics::evaluateMatrix, this);

    workingConditions->updatedReducedField.addListener(&ElectronKinetics::evaluateFieldOperator, this);

    this->eedfType = setup.eedfType;
    this->shapeParameter = setup.shapeParameter;
    this->mixingParameter = setup.numerics.nonLinearRoutines.mixingParameter;
    this->maxEedfRelError = setup.numerics.nonLinearRoutines.maxEedfRelError;
    this->ionizationOperatorType = setup.ionizationOperatorType;
    this->growthModelType = setup.growthModelType;
    this->includeEECollisions = setup.includeEECollisions;
    this->maxPowerBalanceRelError = setup.numerics.maxPowerBalanceRelError;

    // this->plot("Total Elastic Cross Section N2", "Energy (eV)", "Cross Section (m^2)",
    // mixture.grid->getNodes(), mixture.totalCrossSection);

    inelasticMatrix.setZero(grid.cellNumber, grid.cellNumber);

    ionConservativeMatrix.setZero(grid.cellNumber, grid.cellNumber);

    attachmentConservativeMatrix.setZero(grid.cellNumber, grid.cellNumber);

    if (ionizationOperatorType != IonizationOperatorType::conservative &&
        mixture.hasCollisions[static_cast<uint8_t>(CollisionType::ionization)])
    {
        ionizationMatrix.setZero(grid.cellNumber, grid.cellNumber);
    }

    A.setZero(grid.cellNumber);
    B.setZero(grid.cellNumber);

    boltzmannMatrix.setZero();

    // SPARSE INITIALIZATION
    std::vector<Eigen::Triplet<double>> tridiagPattern;
    tridiagPattern.reserve(3 * grid.cellNumber - 2);

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        if (k > 0)
            tridiagPattern.emplace_back(k - 1, k, 0.);

        tridiagPattern.emplace_back(k, k, 0.);

        if (k < grid.cellNumber - 1)
            tridiagPattern.emplace_back(k + 1, k, 0.);
    }

    elasticMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());
    fieldMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());

    if (!mixture.CARGases.empty())
        CARMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());

    /** \todo Could the matrices that are guaranteed to be diagonal just be Eigen::DiagonalMatrix?
     *        That saves a lot of time initializing and makes accessing easier (no 'coeffRef').
     *        Then the pattern-code below can be removed, only a resize needed. Question:
     *        does Eigen optimize matrix addition and multiplication for such matrices? That
     *        is not immediately clear to me.
     *        Same for the other constructor, of course.
     */
    std::vector<Eigen::Triplet<double>> diagPattern;
    diagPattern.reserve(grid.cellNumber);

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        diagPattern.emplace_back(k, k, 0.);
    }

    if (mixture.hasCollisions[static_cast<uint8_t>(CollisionType::attachment)])
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

ElectronKinetics::ElectronKinetics(const json_type &cnf, WorkingConditions *workingConditions)
    : workingConditions(workingConditions), grid(cnf.at("numerics").at("energyGrid")),
      mixture(&grid, cnf, workingConditions), attachmentConservativeMatrix(grid.cellNumber, grid.cellNumber),
      boltzmannMatrix(grid.cellNumber, grid.cellNumber), elasticMatrix(grid.cellNumber, grid.cellNumber),
      fieldMatrix(grid.cellNumber, grid.cellNumber), attachmentMatrix(grid.cellNumber, grid.cellNumber),
      ionSpatialGrowthD(grid.cellNumber, grid.cellNumber), ionSpatialGrowthU(grid.cellNumber, grid.cellNumber),
      fieldMatrixSpatGrowth(grid.cellNumber, grid.cellNumber), fieldMatrixTempGrowth(grid.cellNumber, grid.cellNumber),
      ionTemporalGrowth(grid.cellNumber, grid.cellNumber), g_c(grid.cellNumber), eedf(grid.cellNumber)
{
    grid.updatedMaxEnergy2.addListener(&ElectronKinetics::evaluateMatrix, this);

    workingConditions->updatedReducedField.addListener(&ElectronKinetics::evaluateFieldOperator, this);

    this->eedfType = getEedfType(cnf.at("eedfType"));
    this->shapeParameter = cnf.contains("shapeParameter") ? cnf.at("shapeParameter").get<unsigned>() : 0;
    this->mixingParameter = cnf.at("numerics").at("nonLinearRoutines").at("mixingParameter");
    this->maxEedfRelError = cnf.at("numerics").at("nonLinearRoutines").at("maxEedfRelError");
    this->ionizationOperatorType = getIonizationOperatorType(cnf.at("ionizationOperatorType"));
    this->growthModelType = getGrowthModelType(cnf.at("growthModelType"));
    this->includeEECollisions = cnf.at("includeEECollisions");
    this->maxPowerBalanceRelError = cnf.at("numerics").at("maxPowerBalanceRelError");

    // this->plot("Total Elastic Cross Section N2", "Energy (eV)", "Cross Section (m^2)",
    // mixture.grid->getNodes(), mixture.totalCrossSection);

    inelasticMatrix.setZero(grid.cellNumber, grid.cellNumber);

    ionConservativeMatrix.setZero(grid.cellNumber, grid.cellNumber);

    attachmentConservativeMatrix.setZero(grid.cellNumber, grid.cellNumber);

    if (ionizationOperatorType != IonizationOperatorType::conservative &&
        mixture.hasCollisions[static_cast<uint8_t>(CollisionType::ionization)])
    {
        ionizationMatrix.setZero(grid.cellNumber, grid.cellNumber);
    }

    A.setZero(grid.cellNumber);
    B.setZero(grid.cellNumber);

    boltzmannMatrix.setZero();

    // SPARSE INITIALIZATION
    std::vector<Eigen::Triplet<double>> tridiagPattern;
    tridiagPattern.reserve(3 * grid.cellNumber - 2);

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        if (k > 0)
            tridiagPattern.emplace_back(k - 1, k, 0.);

        tridiagPattern.emplace_back(k, k, 0.);

        if (k < grid.cellNumber - 1)
            tridiagPattern.emplace_back(k + 1, k, 0.);
    }

    elasticMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());
    fieldMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());

    if (!mixture.CARGases.empty())
        CARMatrix.setFromTriplets(tridiagPattern.begin(), tridiagPattern.end());

    std::vector<Eigen::Triplet<double>> diagPattern;
    diagPattern.reserve(grid.cellNumber);

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        diagPattern.emplace_back(k, k, 0.);
    }

    if (mixture.hasCollisions[static_cast<uint8_t>(CollisionType::attachment)])
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

void ElectronKinetics::solve()
{
    if (includeNonConservativeIonization || includeNonConservativeAttachment || includeEECollisions)
    {
        this->mixingDirectSolutions();
    }
    else
    {
        this->invertLinearMatrix();
    }

    if (grid.isSmart)
    {
        double decades = log10(eedf[0]) - log10(eedf[grid.cellNumber - 1]);

        while (decades < grid.minEedfDecay)
        {
            grid.updateMaxEnergy(grid.lastNode() * (1 + grid.updateFactor));

            if (includeNonConservativeIonization || includeNonConservativeAttachment || includeEECollisions)
            {
                this->mixingDirectSolutions();
            }
            else
            {
                this->invertLinearMatrix();
            }

            decades = log10(eedf[0]) - log10(eedf[grid.cellNumber - 1]);
        }

        while (decades > grid.maxEedfDecay)
        {
            grid.updateMaxEnergy(grid.lastNode() / (1 + grid.updateFactor));

            if (includeNonConservativeIonization || includeNonConservativeAttachment || includeEECollisions)
            {
                this->mixingDirectSolutions();
            }
            else
            {
                this->invertLinearMatrix();
            }

            decades = log10(eedf[0]) - log10(eedf[grid.cellNumber - 1]);
        }
    }

    //        for (uint32_t i = 0; i < eedf.size(); ++i) {
    //            printf("%.16e\n", eedf[i]);
    //        }

    evaluatePower(true);

    mixture.evaluateRateCoefficients(eedf);

    evaluateSwarmParameters();

    evaluateFirstAnisotropy();

    obtainedNewEedf.emit(grid, eedf, *workingConditions, power, mixture.gases(), swarmParameters,
                         mixture.rateCoefficients, mixture.rateCoefficientsExtra, firstAnisotropy);

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
    if (!mixture.CARGases.empty())
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
    // Induce normalization condition
    //        matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;
    //        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMatrix = matrix;

    auto begin = std::chrono::high_resolution_clock::now();

    if (!hasSuperelastics)
    {
        eedf.setZero();
        eedf[0] = 1.;

        matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;

        LinAlg::hessenberg(matrix.data(), eedf.data(), grid.cellNumber);
    }
    else
    {
        // TODO: Find a way to distinguish when to use LU and when to use Hessenberg reduction.

        // HESSENBERG WITH PARTIAL PIVOTING
        //            auto *p = new uint32_t[grid.cellNumber];
        //
        //            for (uint32_t i = 0; i < grid.cellNumber; ++i)
        //                p[i] = i;
        //
        //            LinAlg::hessenbergReductionPartialPiv(matrix.data(), &superElasticThresholds[0], p,
        //            grid.cellNumber,
        //                                                  superElasticThresholds.size());
        //
        //            eedf.setZero();
        //            eedf[0] = 1.;
        //
        //            matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;
        //
        //            LinAlg::hessenberg(matrix.data(), eedf.data(), grid.cellNumber);
        //
        //            delete[] p;

        // LU DECOMPOSITION
        Vector b = Vector::Zero(grid.cellNumber);
        b[0] = 1;

        matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;
        eedf = matrix.partialPivLu().solve(b);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "Inverted matrix elapsed time = "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;

    /** \todo It seems that the normaization is superfluous, since the normalization condition
     *        is already part of the system (first row of A, first element of b). One could
     *        decide to change the first equation into eedf[0] = 1 and do the normalization
     *        afterwards. That prevents a fully populated first row of the system matrix
     *        (better sparsity pattern).
     */
    // std::cout << "NORM: " <<  eedf.dot(grid.getCells().cwiseSqrt() * grid.step) << std::endl;
    eedf /= eedf.dot(grid.getCells().cwiseSqrt() * grid.step);
}

void ElectronKinetics::evaluateMatrix()
{
    mixture.evaluateTotalAndElasticCS();

    evaluateElasticOperator();

    evaluateFieldOperator();

    if (!mixture.CARGases.empty())
        evaluateCAROperator();

    evaluateInelasticOperators();

    if (mixture.hasCollisions[static_cast<uint8_t>(CollisionType::ionization)])
        evaluateIonizationOperator();

    if (mixture.hasCollisions[static_cast<uint8_t>(CollisionType::attachment)])
        evaluateAttachmentOperator();

    // Sort and erase duplicates.
    std::sort(superElasticThresholds.begin(), superElasticThresholds.end());
    superElasticThresholds.erase(unique(superElasticThresholds.begin(), superElasticThresholds.end()),
                                 superElasticThresholds.end());
}

void ElectronKinetics::evaluateElasticOperator()
{
    const double Tg = workingConditions->gasTemperature;

    const double factor1 = (Constant::kBeV * Tg / grid.step + 0.5) / grid.step;
    const double factor2 = (Constant::kBeV * Tg / grid.step - 0.5) / grid.step;

    g_c = grid.getNodes().cwiseAbs2().cwiseProduct(mixture.elasticCrossSection) * 2;

    g_c[0] = 0.;
    g_c[g_c.size() - 1] = 0.;

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        elasticMatrix.coeffRef(k, k) = -(g_c[k] * factor1 + g_c[k + 1] * factor2);

        if (k > 0)
            elasticMatrix.coeffRef(k, k - 1) = g_c[k] * factor2;

        if (k < grid.cellNumber - 1)
            elasticMatrix.coeffRef(k, k + 1) = g_c[k + 1] * factor1;
    }
}

void ElectronKinetics::evaluateFieldOperator()
{
    const double EoN = workingConditions->reducedFieldSI;
    const double WoN = workingConditions->reducedExcFreqSI;
    const double me = Constant::electronMass;
    const double e = Constant::electronCharge;

    Vector &cs = mixture.totalCrossSection;

    /** \todo the follosing line produces a NaN for g_E[0] (as expected).
     *        This is later set to 0, so everything is fine. However:
     *        this will crash the code os FPU=exceptions are enabled...
     */
    g_E = ((EoN * EoN / 3) * grid.getNodes()).array() /
          (cs.array() + (me * WoN * WoN / (2 * e)) / (grid.getNodes().cwiseProduct(cs)).array());

    g_E[0] = 0.;
    g_E[g_E.size() - 1] = 0.;

    const double sqStep = grid.step * grid.step;

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        fieldMatrix.coeffRef(k, k) = -(g_E[k] + g_E[k + 1]) / sqStep;

        if (k > 0)
            fieldMatrix.coeffRef(k, k - 1) = g_E[k] / sqStep;

        if (k < grid.cellNumber - 1)
            fieldMatrix.coeffRef(k, k + 1) = g_E[k + 1] / sqStep;
    }
}

void ElectronKinetics::evaluateCAROperator()
{
    const double Tg = workingConditions->gasTemperature;

    const double factor1 = (Constant::kBeV * Tg / grid.step + 0.5) / grid.step;
    const double factor2 = (Constant::kBeV * Tg / grid.step - 0.5) / grid.step;

    double sigma0B = 0.;

    for (auto &gas : mixture.CARGases)
    {
        sigma0B += gas->fraction * gas->electricQuadrupoleMoment * gas->rotationalConstant;
    }

    sigma0B *= 8. * Constant::pi / (15. * Constant::electronCharge);
    g_CAR = grid.getNodes() * (4. * sigma0B);

    // Boundary conditions
    g_CAR[0] = 0.;
    g_CAR[grid.cellNumber] = 0.;

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        CARMatrix.coeffRef(k, k) = -(g_CAR[k] * factor1 + g_CAR[k + 1] * factor2);

        if (k > 0)
            CARMatrix.coeffRef(k, k - 1) = g_CAR[k] * factor2;

        if (k < grid.cellNumber - 1)
            CARMatrix.coeffRef(k, k + 1) = g_CAR[k + 1] * factor1;
    }
}

void ElectronKinetics::evaluateInelasticOperators()
{
    const uint32_t cellNumber = grid.cellNumber;

    inelasticMatrix.setZero();

    for (const auto &gas : mixture.gases())
    {
        for (auto vecIndex = static_cast<uint8_t>(CollisionType::excitation);
             vecIndex <= static_cast<uint8_t>(CollisionType::rotational); ++vecIndex)
        {

            for (const auto &collision : gas->collisions[vecIndex])
            {
                const double threshold = collision->crossSection->threshold;

                if (threshold < grid.step || threshold > grid.getNodes()[grid.cellNumber])
                    continue;

                const double targetDensity = collision->getTarget()->density;

                if (targetDensity != 0)
                {
                    const auto numThreshold = static_cast<uint32_t>(std::floor(threshold / grid.step));

                    Vector cellCrossSection(cellNumber);

                    for (uint32_t i = 0; i < cellNumber; ++i)
                        cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

                    for (uint32_t k = 0; k < cellNumber; ++k)
                    {
                        if (k < cellNumber - numThreshold)
                            inelasticMatrix(k, k + numThreshold) +=
                                targetDensity * grid.getCells()[k + numThreshold] * cellCrossSection[k + numThreshold];

                        inelasticMatrix(k, k) -= targetDensity * grid.getCells()[k] * cellCrossSection[k];
                    }

                    if (collision->isReverse)
                    {
                        const double swRatio = collision->getTarget()->statisticalWeight /
                                               collision->m_rhsHeavyStates[0]->statisticalWeight;
                        const double productDensity = collision->m_rhsHeavyStates[0]->density;

                        if (productDensity == 0)
                            continue;

                        if (numThreshold > 1)
                            superElasticThresholds.emplace_back(numThreshold);

                        if (numThreshold != 1)
                            hasSuperelastics = true;

                        for (uint32_t k = 0; k < cellNumber; ++k)
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

    for (const auto &gas : mixture.gases())
    {
        for (const auto &collision : gas->collisions[static_cast<uint8_t>(CollisionType::ionization)])
        {
            const double threshold = collision->crossSection->threshold;

            if (threshold > grid.getNode(grid.cellNumber))
                continue;

            hasValidCollisions = true;

            const double density = collision->getTarget()->density;
            const auto numThreshold = static_cast<uint32_t>(std::floor(threshold / grid.step));

            Vector cellCrossSection(grid.cellNumber);

            for (uint32_t i = 0; i < grid.cellNumber; ++i)
                cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

            switch (ionizationOperatorType)
            {
            case IonizationOperatorType::conservative:
                break;

            case IonizationOperatorType::oneTakesAll:
                for (uint32_t k = 0; k < grid.cellNumber; ++k)
                {
                    if (k < grid.cellNumber - numThreshold)
                        ionizationMatrix(k, k + numThreshold) +=
                            density * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                    const double term = density * grid.getCell(k) * cellCrossSection(k);

                    ionizationMatrix(k, k) -= term;
                    ionizationMatrix(0, k) += term;
                }
                break;

            case IonizationOperatorType::equalSharing:
                for (uint32_t k = 0; k < grid.cellNumber; ++k)
                {
                    ionizationMatrix(k, k) -= density * grid.getCell(k) * cellCrossSection[k];

                    if (k < (grid.cellNumber - numThreshold) / 2)
                    {
                        const uint32_t i = 2 * (k + 1) + numThreshold - 1;

                        ionizationMatrix(k, i) += 4 * density * grid.getCell(i) * cellCrossSection(i);
                    }
                }
                break;
            case IonizationOperatorType::sdcs:
                double W = gas->OPBParameter;

                if (W < 0)
                    W = threshold;

                for (uint32_t k = 0; k < grid.cellNumber; ++k)
                {
                    const uint32_t end = std::min(2 * (k + 1) + numThreshold, grid.cellNumber);

                    if (k > numThreshold)
                    {
                        const uint32_t half = (k + 1 - numThreshold) / 2;
                        const double numerator = 1 / std::atan((grid.getCell(k) - threshold) / (2 * W));
                        double sum = 0.;

                        for (uint32_t i = 0; i < half; ++i)
                            sum += numerator / (W + grid.getCell(i) * grid.getCell(i) / W);

                        ionizationMatrix(k, k) -= density * grid.step * grid.getCell(k) * cellCrossSection[k] * sum;
                    }

                    /** \todo If k + numThreshold + 1 < grid.cellNumber, the term is ignored.
                     *        Document (in the document, not necessarily here) what are the
                     *        consequences of that.
                     */
                    if (k + numThreshold + 1 < grid.cellNumber)
                    {
                        for (uint32_t i = k + numThreshold + 1; i < end; ++i)
                        {
                            ionizationMatrix(k, i) += density * grid.step * grid.getCell(i) * cellCrossSection[i] /
                                                      (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                       (W + std::pow(grid.getCell(i - k - numThreshold - 1), 2) / W));
                        }
                    }

                    /** \todo The following comment needs to be sorted out (possible index errors).
                     *  \todo Document what is done here (algorithm) and how it is implemented.
                     */

                    // This last section might need some adjustments because of indexing
                    // differences between Matlab and C++ (since indexes are multiplied here).

                    for (uint32_t i = 2 * (k + 1) + numThreshold - 1; i < grid.cellNumber; ++i)
                    {
                        ionizationMatrix(k, i) += density * grid.step * grid.getCell(i) * cellCrossSection[i] /
                                                  (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                                   (W + std::pow(grid.getCell(k), 2) / W));
                    }
                }
                break;
            }

            // Evaluation of the conservative ionization operator

            if (numThreshold == 0)
                continue;

            for (uint32_t k = 0; k < grid.cellNumber; ++k)
            {
                if (k < grid.cellNumber - numThreshold)
                    ionConservativeMatrix(k, k + numThreshold) +=
                        density * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                ionConservativeMatrix(k, k) -= density * grid.getCell(k) * cellCrossSection[k];
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

    const uint32_t cellNumber = grid.cellNumber;

    for (const auto &gas : mixture.gases())
    {
        for (const auto &collision : gas->collisions[static_cast<uint8_t>(CollisionType::attachment)])
        {
            const double threshold = collision->crossSection->threshold;

            if (threshold > grid.getNode(cellNumber))
                continue;

            /** \this should definitely not be in this (double) loop. Is this a constructor task?
             *        Can this just be replaced with 'gas->collisions[CollisionType::attachment].size()'
             *        in (other) places where this is now used?
             */
            includeNonConservativeAttachment = true;

            const auto numThreshold = static_cast<uint32_t>(std::floor(threshold / grid.step));

            /** \todo Eliminate the cellCrossSection vector? This can be calculated on the fly
             *        in the two places where it is needed (one if the merger below can be done).
             */
            Vector cellCrossSection(cellNumber);

            const double targetDensity = collision->getTarget()->density;

            /// \todo Merge with the subsequent k-loop.
            for (uint32_t i = 0; i < cellNumber; ++i)
                cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

            for (uint32_t k = 0; k < cellNumber; ++k)
                attachmentMatrix.coeffRef(k, k) -= targetDensity * grid.getCell(k) * cellCrossSection[k];

            if (numThreshold == 0)
                continue;

            /** Can we also merge with this loop?
             *  It does not seem problematic to have an 'if (numThreshold)' in the loop,
             *  given the amount of work done.
             */
            for (uint32_t k = 0; k < cellNumber; ++k)
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

    const uint32_t numCells = grid.cellNumber;

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

                if (((eedfOld - eedf).cwiseAbs().array() / eedfOld.array()).maxCoeff() < maxEedfRelError)
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
    const double e = Constant::electronCharge, m = Constant::electronMass, EoN = workingConditions->reducedFieldSI;

    Vector cellTotalCrossSection(grid.cellNumber);

    for (uint32_t i = 0; i < grid.cellNumber; ++i)
        cellTotalCrossSection[i] = .5 * (mixture.totalCrossSection[i] + mixture.totalCrossSection[i + 1]);

    if (!mixture.CARGases.empty())
    {
        boltzmannMatrix =
            1.e20 * (elasticMatrix + fieldMatrix + CARMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
    }
    else
    {
        boltzmannMatrix = 1.e20 * (elasticMatrix + fieldMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
    }

    Vector baseDiag(grid.cellNumber), baseSubDiag(grid.cellNumber), baseSupDiag(grid.cellNumber);

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        baseDiag[k] = boltzmannMatrix(k, k);

        if (k > 0)
            baseSubDiag[k] = boltzmannMatrix(k, k - 1);

        if (k < grid.cellNumber - 1)
            baseSupDiag[k] = boltzmannMatrix(k, k + 1);
    }

    if (includeEECollisions)
    {
        A = alphaEE / grid.step * (BAee.transpose() * eedf);
        B = alphaEE / grid.step * (BAee * eedf);
    }

    Vector integrandCI = (sqrt(2. * e / m) * grid.step) * Vector::Ones(grid.cellNumber).transpose() *
                         (ionizationMatrix + attachmentMatrix);

    double CIEffNew = eedf.dot(integrandCI);
    double CIEffOld = CIEffNew / 3;

    CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

    // diffusion and mobility components of the spatial growth terms
    // can be removed since this is already done in the directMixing function
    ionSpatialGrowthD.setZero();
    ionSpatialGrowthU.setZero();

    Vector tempVector = grid.getCells().array() / (3. * cellTotalCrossSection).array();

    Vector D0 = tempVector, U0sup(grid.cellNumber), U0inf(grid.cellNumber);

    U0sup[0] = 0.;
    U0inf[grid.cellNumber - 1] = 0.;

    for (uint32_t j = 0; j < grid.cellNumber; ++j)
    {
        if (j != 0)
            U0sup[j] = EoN / (2. * grid.step) * tempVector[j - 1];

        if (j != grid.cellNumber - 1)
            U0inf[j] = -EoN / (2. * grid.step) * tempVector[j + 1];
    }

    Vector U0 = U0sup + U0inf;

    double ND = sqrt(2 * e / m) * grid.step * D0.dot(eedf), muE = -sqrt(2 * e / m) * grid.step * U0.dot(eedf);

    double alphaRedEffNew, alphaRedEffOld = 0.;

    if (muE * muE - 4 * CIEffNew * ND < 0.)
    {
        alphaRedEffNew = CIEffNew / muE;
    }
    else
    {
        alphaRedEffNew = (muE - sqrt(muE * muE - 4 * CIEffNew * ND)) / (2 * ND);
    }

    uint32_t iter = 0;
    bool hasConverged = false;

    const Vector g_fieldSpatialBase = (EoN / 6) * grid.getNodes().array() / mixture.totalCrossSection.array();

    while (!hasConverged)
    {
        g_fieldSpatialGrowth = alphaRedEffNew * g_fieldSpatialBase;
        g_fieldSpatialGrowth[0] = 0.;
        g_fieldSpatialGrowth[grid.cellNumber] = 0.;

        for (uint32_t k = 0; k < grid.cellNumber; ++k)
        {
            fieldMatrixSpatGrowth.coeffRef(k, k) = (g_fieldSpatialGrowth[k + 1] - g_fieldSpatialGrowth[k]) / grid.step;

            if (k > 0)
                fieldMatrixSpatGrowth.coeffRef(k, k - 1) = -g_fieldSpatialGrowth[k] / grid.step;

            if (k < grid.cellNumber - 1)
                fieldMatrixSpatGrowth.coeffRef(k, k + 1) = g_fieldSpatialGrowth[k + 1] / grid.step;
        }

        for (uint32_t k = 0; k < grid.cellNumber; ++k)
        {
            ionSpatialGrowthD.coeffRef(k, k) = alphaRedEffNew * alphaRedEffNew * D0[k];

            if (k > 0)
                ionSpatialGrowthU.coeffRef(k, k - 1) = alphaRedEffNew * U0inf[k - 1];

            if (k < grid.cellNumber - 1)
                ionSpatialGrowthU.coeffRef(k, k + 1) = alphaRedEffNew * U0sup[k + 1];
        }

        for (uint32_t k = 0; k < grid.cellNumber; ++k)
        {
            boltzmannMatrix(k, k) = baseDiag[k] + 1.e20 * (fieldMatrixSpatGrowth.coeff(k, k) +
                                                           ionSpatialGrowthD.coeff(k, k) - (A[k] + B[k]));

            if (k > 0)
                boltzmannMatrix(k, k - 1) = baseSubDiag[k] + 1.e20 * (fieldMatrixSpatGrowth.coeff(k, k - 1) +
                                                                      ionSpatialGrowthU.coeff(k, k - 1) + A[k - 1]);

            if (k < grid.cellNumber - 1)
                boltzmannMatrix(k, k + 1) = baseSupDiag[k] + 1.e20 * (fieldMatrixSpatGrowth.coeff(k, k + 1) +
                                                                      ionSpatialGrowthU.coeff(k, k + 1) + B[k + 1]);
        }

        Vector eedfNew = eedf;

        invertMatrix(boltzmannMatrix);

        CIEffOld = CIEffNew;
        CIEffNew = eedf.dot(integrandCI);

        CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

        ND = sqrt(2 * e / m) * grid.step * D0.dot(eedf);
        muE = -sqrt(2 * e / m) * grid.step * U0.dot(eedf);

        const double discriminant = muE * muE - 4 * CIEffNew * ND;

        alphaRedEffOld = alphaRedEffNew;
        alphaRedEffNew = (discriminant < 0) ? CIEffNew / muE : (muE - sqrt(discriminant)) / (2 * ND);

        alphaRedEffNew = mixingParameter * alphaRedEffNew + (1 - mixingParameter) * alphaRedEffOld;

        if (((alphaRedEffNew == 0 || abs(alphaRedEffNew - alphaRedEffOld) / alphaRedEffOld < 1.e-10) &&
             ((eedf - eedfNew).cwiseAbs().array() / eedf.array()).maxCoeff() < maxEedfRelError) ||
            iter > 150)
        {
            hasConverged = true;

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
    const double e = Constant::electronCharge, m = Constant::electronMass, EoN = workingConditions->reducedFieldSI,
                 WoN = workingConditions->reducedExcFreqSI;

    if (!mixture.CARGases.empty())
    {
        boltzmannMatrix = 1.e20 * (elasticMatrix + CARMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
    }
    else
    {
        boltzmannMatrix = 1.e20 * (elasticMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
    }

    Vector baseDiag(grid.cellNumber), baseSubDiag(grid.cellNumber), baseSupDiag(grid.cellNumber);

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        baseDiag[k] = boltzmannMatrix(k, k);

        if (k > 0)
            baseSubDiag[k] = boltzmannMatrix(k, k - 1);

        if (k < grid.cellNumber - 1)
            baseSupDiag[k] = boltzmannMatrix(k, k + 1);
    }

    if (includeEECollisions)
    {
        A = alphaEE / grid.step * (BAee.transpose() * eedf);
        B = alphaEE / grid.step * (BAee * eedf);
    }

    Vector integrandCI = (sqrt(2. * e / m) * grid.step) * Vector::Ones(grid.cellNumber).transpose() *
                         (ionizationMatrix + attachmentMatrix);
    ;

    double CIEffNew = eedf.dot(integrandCI);
    double CIEffOld = CIEffNew / 3.;

    CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

    Vector totalCSI(grid.cellNumber + 1), eedfNew(grid.cellNumber);

    bool hasConverged = false;
    uint32_t iter = 0;

    while (!hasConverged)
    {
        Log<Message>::Notify("Iteration ", iter);

        totalCSI[0] = mixture.totalCrossSection[0];

        const long double growthFactor = CIEffNew * sqrt(m / (2 * e));

        for (uint32_t i = 1; i <= grid.cellNumber; ++i)
        {
            totalCSI[i] = mixture.totalCrossSection[i] + growthFactor / sqrt(i * grid.step);
        }

        g_fieldTemporalGrowth =
            ((EoN * EoN / 3) * grid.getNodes()).array() /
            (totalCSI.array() + (m * WoN * WoN / (2 * e)) / (grid.getNodes().cwiseProduct(totalCSI)).array());
        g_fieldTemporalGrowth[0] = 0.;
        g_fieldTemporalGrowth[grid.cellNumber] = 0.;

        const double sqrStep = grid.step * grid.step;

        for (uint32_t k = 0; k < grid.cellNumber; ++k)
        {
            fieldMatrixTempGrowth.coeffRef(k, k) = -(g_fieldTemporalGrowth[k] + g_fieldTemporalGrowth[k + 1]) / sqrStep;

            if (k > 0)
                fieldMatrixTempGrowth.coeffRef(k, k - 1) = g_fieldTemporalGrowth[k] / sqrStep;

            if (k < grid.cellNumber - 1)
                fieldMatrixTempGrowth.coeffRef(k, k + 1) = g_fieldTemporalGrowth[k + 1] / sqrStep;

            ionTemporalGrowth.coeffRef(k, k) = -growthFactor * sqrt(grid.getCell(k));
        }

        for (uint32_t k = 0; k < grid.cellNumber; ++k)
        {
            boltzmannMatrix(k, k) = baseDiag[k] + 1.e20 * (fieldMatrixTempGrowth.coeff(k, k) +
                                                           ionTemporalGrowth.coeff(k, k) - (A[k] + B[k]));

            if (k > 0)
                boltzmannMatrix(k, k - 1) = baseSubDiag[k] + 1.e20 * (fieldMatrixTempGrowth.coeff(k, k - 1) + A[k - 1]);

            if (k < grid.cellNumber - 1)
                boltzmannMatrix(k, k + 1) = baseSupDiag[k] + 1.e20 * (fieldMatrixTempGrowth.coeff(k, k + 1) + B[k + 1]);
        }

        eedfNew = eedf;

        invertMatrix(boltzmannMatrix);

        CIEffOld = CIEffNew;
        CIEffNew = eedf.dot(integrandCI);
        CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

        if (((CIEffNew == 0 || abs(CIEffNew - CIEffOld) / CIEffOld < 1.e10) &&
             ((eedf - eedfNew).cwiseAbs().array() / eedf.array()).maxCoeff() < maxEedfRelError) ||
            iter > 150)
        {
            hasConverged = true;

            if (iter > 150 && !includeEECollisions)
                Log<Message>::Warning("Iterative temporal growth scheme did not converge.");
        }

        ++iter;
    }

    std::cerr << "Temporal growth routine converged in: " << iter << " iterations.\n";

    CIEff = CIEffOld;
}

void ElectronKinetics::solveEEColl()
{
    const double e = Constant::electronCharge, e0 = Constant::vacuumPermittivity,
                 ne = workingConditions->electronDensity, n0 = workingConditions->gasDensity;

    // Splitting all possible options for the best performance.

    if (includeNonConservativeIonization || includeNonConservativeAttachment)
    {
        if (growthModelType == GrowthModelType::spatial)
        {
            if (mixture.CARGases.empty())
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
            if (mixture.CARGases.empty())
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
        if (mixture.CARGases.empty())
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

    Vector baseDiag(grid.cellNumber), baseSubDiag(grid.cellNumber), baseSupDiag(grid.cellNumber);

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        baseDiag[k] = boltzmannMatrix(k, k);

        if (k > 0)
            baseSubDiag[k] = boltzmannMatrix(k, k - 1);

        if (k < grid.cellNumber - 1)
            baseSupDiag[k] = boltzmannMatrix(k, k + 1);
    }

    BAee.setZero(grid.cellNumber, grid.cellNumber);

    const Vector cellsThreeOverTwo = grid.getCells().cwiseProduct(grid.getCells().cwiseSqrt()),
                 energyArray = -(grid.step / 2.) * grid.getCells().cwiseSqrt() + (2. / 3.) * cellsThreeOverTwo;

    for (uint32_t j = 0; j < grid.cellNumber - 1; ++j)
    {

        for (uint32_t i = 1; i <= j; ++i)
            BAee(i, j) = energyArray[i];

        const double value = 2. / 3. * std::pow(grid.getNode(j + 1), 1.5);

        for (uint32_t i = j + 1; i < grid.cellNumber; ++i)
            BAee(i, j) = value;
    }

    // detailed balance condition

    for (uint32_t j = 0; j < grid.cellNumber - 1; ++j)
    {
        for (uint32_t i = 1; i < grid.cellNumber; ++i)
        {
            BAee(i, j) = sqrt(BAee(i, j) * BAee(j + 1, i - 1));
        }
    }

    double meanEnergy = grid.step * cellsThreeOverTwo.dot(eedf), Te = 2. / 3. * meanEnergy,
           logC = std::log(12 * Constant::pi * std::pow(e0 * Te / e, 1.5) / std::sqrt(ne)),
           alpha = (ne / n0) * (e * e / (8 * Constant::pi * e0 * e0)) * logC;

    double ratioNew = 0.;
    Vector eedfNew = eedf;

    bool hasConverged = false;
    uint32_t iter = 0;

    //        Vector MeeDiag(grid.cellNumber), MeeSub(grid.cellNumber), MeeSup(grid.cellNumber);

    // In this implementation we completely skip the Mee matrix, saving both memory and time.

    while (!hasConverged)
    {
        A = (alpha / grid.step) * (BAee.transpose() * eedf);
        B = (alpha / grid.step) * (BAee * eedf);

        for (uint32_t k = 0; k < grid.cellNumber; ++k)
        {
            boltzmannMatrix(k, k) = baseDiag[k] - 1.e20 * (A[k] + B[k]);

            if (k > 0)
                boltzmannMatrix(k, k - 1) = baseSubDiag[k] + 1.e20 * A[k - 1];

            if (k < grid.cellNumber - 1)
                boltzmannMatrix(k, k + 1) = baseSupDiag[k] + 1.e20 * B[k + 1];
        }

        invertMatrix(boltzmannMatrix);

        evaluatePower(false);

        const double ratio = std::abs(power.electronElectron / power.reference);

        if (((eedf - eedfNew).cwiseAbs().array() / eedf.array()).maxCoeff() < maxEedfRelError)
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

            Vector eedfOld = eedfNew;
            eedfNew = eedf;

            eedf = eedfNew - (ratioNew / (ratioNew - ratioOld)) * (eedfNew - eedfOld);

            for (uint32_t i = 0; i < grid.cellNumber; ++i)
            {
                if (eedf[i] < 0)
                    eedf[i] = abs(eedf[i]);
            }
        }

        meanEnergy = grid.step * cellsThreeOverTwo.dot(eedf);
        Te = 2. / 3. * meanEnergy;
        logC = std::log(12 * Constant::pi * std::pow(e0 * Te / e, 1.5) / std::sqrt(ne));
        alpha = (ne / n0) * (e * e / (8 * Constant::pi * e0 * e0)) * logC;

        iter++;
    }

    alphaEE = alpha;

    std::cerr << "e-e routine converged in: " << iter << " iterations.\n";
}

void ElectronKinetics::evaluatePower(bool isFinalSolution)
{
    const double factor = sqrt(2. * Constant::electronCharge / Constant::electronMass),
                 kTg = Constant::kBeV * workingConditions->gasTemperature,
                 auxHigh = kTg + grid.step * .5, // aux1
        auxLow = kTg - grid.step * .5;           // aux2

    power = Power();

    double elasticNet = 0., elasticGain = 0.;

    for (uint32_t k = 0; k < grid.cellNumber; ++k)
    {
        elasticNet += eedf[k] * (g_c[k + 1] * auxLow - g_c[k] * auxHigh);
        elasticGain += eedf[k] * (g_c[k + 1] - g_c[k]);
    }

    power.elasticNet = factor * elasticNet;
    power.elasticGain = factor * kTg * elasticGain;
    power.elasticLoss = power.elasticNet - power.elasticGain;

    if (!mixture.CARGases.empty())
    {
        double carNet = 0., carGain = 0.;

        for (uint32_t k = 0; k < grid.cellNumber - 1; ++k)
        {
            carNet += eedf[k] * (g_CAR[k + 1] * auxLow - g_CAR[k] * auxHigh);
            carGain += eedf[k] * (g_CAR[k + 1] - g_CAR[k]);
        }

        power.carNet = factor * carNet;
        power.carGain = factor * kTg * carGain;
        power.carLoss = power.carNet - power.carGain;
    }

    if (includeNonConservativeIonization || includeNonConservativeAttachment)
    {
        if (growthModelType == GrowthModelType::temporal)
        {
            double field = 0., growthModel = 0.;

            for (uint32_t k = 0; k < grid.cellNumber; ++k)
            {
                field += eedf[k] * (g_fieldTemporalGrowth[k + 1] - g_fieldTemporalGrowth[k]);
                growthModel += eedf[k] * grid.getCell(k) * sqrt(grid.getCell(k));
            }

            power.field = factor * field;
            power.eDensGrowth = -CIEff * grid.step * growthModel;
        }
        else if (growthModelType == GrowthModelType::spatial)
        {
            double field = 0., correction = 0., powerDiffusion = 0., powerMobility = 0.;
            Vector cellCrossSection(grid.cellNumber);

            for (uint32_t k = 0; k < grid.cellNumber; ++k)
            {
                field += eedf[k] * (g_E[k + 1] - g_E[k]);
                correction -= eedf[k] * (g_fieldSpatialGrowth[k + 1] + g_fieldSpatialGrowth[k]);

                // Diffusion and Mobility contributions
                cellCrossSection[k] = .5 * (mixture.totalCrossSection[k] + mixture.totalCrossSection[k + 1]);
                powerDiffusion += grid.getCell(k) * grid.getCell(k) * eedf[k] / cellCrossSection[k];

                if (k > 0 && k < grid.cellNumber - 1)
                {
                    powerMobility +=
                        grid.getCell(k) * grid.getCell(k) * (eedf[k + 1] - eedf[k - 1]) / cellCrossSection[k];
                }
            }

            power.field = factor * (field + grid.step * correction);
            power.eDensGrowth = alphaRedEff * alphaRedEff * factor * grid.step / 3. * powerDiffusion +
                                factor * alphaRedEff * (workingConditions->reducedFieldSI / 6.) *
                                    (grid.getCell(0) * grid.getCell(0) * eedf[1] / cellCrossSection[0] -
                                     grid.getCell(grid.cellNumber - 1) * grid.getCell(grid.cellNumber - 1) *
                                         eedf[grid.cellNumber - 2] / cellCrossSection[grid.cellNumber - 1] +
                                     powerMobility);
        }
    }
    else
    {
        double field = 0;

        for (uint32_t k = 0; k < grid.cellNumber; ++k)
        {
            field += eedf[k] * (g_E[k + 1] - g_E[k]);
        }

        power.field = factor * field;
    }

    if (includeEECollisions)
    {
        power.electronElectron = (-factor * grid.step * grid.step) * (A - B).dot(eedf);
    }

    // Evaluate power absorbed per electron at unit gas density due to in- and superelastic collisions.
    for (auto &gas : mixture.gases())
    {
        gas->evaluatePower(ionizationOperatorType, eedf);
        power += gas->getPower();
    }

    power.inelastic =
        power.excitationIne + power.vibrationalIne + power.rotationalIne + power.ionizationIne + power.attachmentIne;
    power.superelastic = power.excitationSup + power.vibrationalSup + power.rotationalSup;

    double totalGain = 0., totalLoss = 0.;

    double powerValues[13]{power.field,           power.elasticGain,   power.elasticLoss,   power.carGain,
                           power.carLoss,         power.excitationSup, power.excitationIne, power.vibrationalSup,
                           power.vibrationalIne,  power.rotationalSup, power.rotationalIne, power.eDensGrowth,
                           power.electronElectron};

    for (double value : powerValues)
    {
        if (value > 0)
            totalGain += value;
        else
            totalLoss += value;
    }

    power.balance = power.field + power.elasticNet + power.carNet + power.inelastic + power.superelastic +
                    power.eDensGrowth + power.electronElectron;
    power.relativeBalance = abs(power.balance) / totalGain;
    power.reference = totalGain;

    if (isFinalSolution && power.relativeBalance > maxPowerBalanceRelError)
        Log<PowerBalanceError>::Warning(maxPowerBalanceRelError);
}

void ElectronKinetics::evaluateSwarmParameters()
{
    const double me = Constant::electronMass, e = Constant::electronCharge;

    const uint32_t n = grid.cellNumber;

    const bool nonConservative = (includeNonConservativeIonization || includeNonConservativeAttachment);

    Vector tCS(mixture.totalCrossSection);

    if (growthModelType == GrowthModelType::temporal && nonConservative)
    {

        tCS.tail(grid.cellNumber).array() +=
            CIEff * std::sqrt(me / (2 * e)) / grid.getNodes().tail(n).cwiseSqrt().array();
    }

    swarmParameters.redDiffCoeff = 2. / 3. * std::sqrt(2. * e / me) * grid.step *
                                   grid.getCells().cwiseProduct(eedf).cwiseQuotient(tCS.head(n) + tCS.tail(n)).sum();

    swarmParameters.redMobCoeff = -std::sqrt(2. * e / me) / 3. *
                                  grid.getNodes()
                                      .segment(1, n - 1)
                                      .cwiseProduct(eedf.tail(n - 1) - eedf.head(n - 1))
                                      .cwiseQuotient(tCS.segment(1, n - 1))
                                      .sum();

    if (growthModelType == GrowthModelType::spatial && nonConservative)
    {
        swarmParameters.driftVelocity = -swarmParameters.redDiffCoeff * alphaRedEff +
                                        swarmParameters.redMobCoeff * workingConditions->reducedFieldSI;
    }
    else
    {
        swarmParameters.driftVelocity = swarmParameters.redMobCoeff * workingConditions->reducedFieldSI;
    }

    double totalIonRateCoeff = 0., totalAttRateCoeff = 0.;

    for (const auto &gas : mixture.gases())
    {
        for (const auto &collision : gas->collisions[static_cast<uint8_t>(CollisionType::ionization)])
        {
            totalIonRateCoeff += collision->getTarget()->density * collision->ineRateCoeff;
        }

        for (const auto &collision : gas->collisions[static_cast<uint8_t>(CollisionType::attachment)])
        {
            totalAttRateCoeff += collision->getTarget()->density * collision->ineRateCoeff;
        }
    }

    swarmParameters.redTownsendCoeff = totalIonRateCoeff / swarmParameters.driftVelocity;
    swarmParameters.redAttCoeff = totalAttRateCoeff / swarmParameters.driftVelocity;

    swarmParameters.meanEnergy = grid.step * (grid.getCells().array().pow(1.5) * eedf.array()).sum();

    swarmParameters.characEnergy = swarmParameters.redDiffCoeff / swarmParameters.redMobCoeff;

    swarmParameters.Te = 2. / 3. * swarmParameters.meanEnergy;

    // TODO: is this correct? (simulations after the first will have a different value for Te).
    workingConditions->updateElectronTemperature(swarmParameters.Te);
}

void ElectronKinetics::evaluateFirstAnisotropy()
{
    firstAnisotropy.setZero(grid.cellNumber);

    const double e = Constant::electronCharge, me = Constant::electronMass, EoN = workingConditions->reducedFieldSI,
                 WoN = workingConditions->reducedExcFreqSI;

    const uint32_t n = grid.cellNumber;

    firstAnisotropy[0] = (eedf[1] - eedf[0]) / grid.step;
    firstAnisotropy[n - 1] = (eedf[n - 1] - eedf[n - 2]) / grid.step;

    firstAnisotropy.segment(1, n - 2) = (eedf.segment(2, n - 2) - eedf.segment(0, n - 2)) / (2 * grid.step);

    Vector cellCrossSection = (mixture.totalCrossSection.segment(0, n) + mixture.totalCrossSection.segment(1, n)) / 2.;

    if (includeNonConservativeIonization || includeNonConservativeAttachment)
    {
        if (growthModelType == GrowthModelType::temporal)
        {
            cellCrossSection =
                cellCrossSection.array() + CIEff / (std::sqrt(e * e / me) * grid.getCells().cwiseSqrt()).array();

            if (WoN == 0)
            {
                firstAnisotropy = -EoN * firstAnisotropy.cwiseQuotient(cellCrossSection);
            }
            else
            {
                firstAnisotropy = -EoN * sqrt(2.) * firstAnisotropy.array() /
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
            -EoN * sqrt(2.) * firstAnisotropy.array() /
            (cellCrossSection.array() + WoN * WoN * me / (2. * e * grid.getCells().array() * cellCrossSection.array()));
    }
}
} // namespace loki

// Code to print a matrix
//            for (uint32_t i = 0; i < grid.cellNumber; ++i) {
//                for (uint32_t j = 0; j < grid.cellNumber; ++j) {
//                    printf("%.16e", matrix(i, j));
//
//                    if (j < grid.cellNumber - 1) {
//                        printf("\t");
//                    }
//                }
//
//                printf("\n");
//            }
