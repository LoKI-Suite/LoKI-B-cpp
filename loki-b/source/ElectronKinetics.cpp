//
// Created by daan on 13-5-19.
//

#include "ElectronKinetics.h"
#include <chrono>
#include <cmath>

// TODO: With the introduction of the job system, the only matrix that we have to be able
//  to update separately is the fieldMatrix. Therefore we can simply store the other matrices
//  as one matrix. This saves a large amount of memory, and, in the case of large matrices,
//  also quite a bit of running time.
//  This is not completely true however, since we need a different combination of matrices
//  in different parts of the simulation.
//  A better solution might be to write a tridiagonal matrix class that stores the elements
//  in three separate vectors. It would be best to do this in an Eigen compliant way, such
//  that these matrices can simply be added to dense matrices (the + and [] operators are 
//  the only operators to overload).

namespace loki {
    ElectronKinetics::ElectronKinetics(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions)
            : workingConditions(workingConditions), grid(setup.numerics.energyGrid), mixture(&grid),
              g_c(grid.cellNumber), eedf(grid.cellNumber), elasticMatrix(grid.cellNumber, grid.cellNumber),
              continuousMatrix(grid.cellNumber, grid.cellNumber), attachmentMatrix(grid.cellNumber, grid.cellNumber),
              attachmentConservativeMatrix(grid.cellNumber, grid.cellNumber) {

        mixture.initialize(setup, workingConditions);

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

        elasticMatrix.setZero(grid.cellNumber, grid.cellNumber);
        fieldMatrix.setZero(grid.cellNumber, grid.cellNumber);
        inelasticMatrix.setZero(grid.cellNumber, grid.cellNumber);

        if (!mixture.CARGasses.empty())
            CARMatrix.setZero(grid.cellNumber, grid.cellNumber);

        ionConservativeMatrix.setZero(grid.cellNumber, grid.cellNumber);

        if (ionizationOperatorType != IonizationOperatorType::conservative)
            ionizationMatrix.setZero(grid.cellNumber, grid.cellNumber);

        attachmentMatrix.setZero(grid.cellNumber, grid.cellNumber);

        this->evaluateMatrix();
    }

    void ElectronKinetics::solve() {
        // TODO: add attachment
        if (includeNonConservativeIonization) {
            this->mixingDirectSolutions();
        } else {
            this->invertLinearMatrix();
        }

        for (uint32_t i = 0; i < eedf.size(); ++i) {
            printf("%.16e\n", eedf[i]);
        }

//        this->plot("Eedf", "Energy (eV)", "Eedf (Au)", grid.getCells(), eedf);
    }

    void ElectronKinetics::invertLinearMatrix() {
        Matrix boltzmannMatrix; // = (elasticMatrix + fieldMatrix + CARMatrix + inelasticMatrix + ionConservativeMatrix).array() * 1.e20;

        if (!mixture.CARGasses.empty()) {
            boltzmannMatrix =
                    (elasticMatrix + fieldMatrix + CARMatrix + inelasticMatrix + ionConservativeMatrix +
                     attachmentConservativeMatrix).array() * 1.e20;
        } else {
            boltzmannMatrix = (elasticMatrix + fieldMatrix + inelasticMatrix + ionConservativeMatrix +
                               attachmentConservativeMatrix).array() * 1.e20;
        }

        invertMatrix(boltzmannMatrix);
    }

    void ElectronKinetics::invertMatrix(Matrix &matrix) {

        // Induce normalization condition
//        matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;
//        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMatrix = matrix;

        auto begin = std::chrono::high_resolution_clock::now();

        if (!hasSuperelastics && !includeEECollisions) {
            eedf.setZero();
            eedf[0] = 1.;

            matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;

            LinAlg::hessenberg(matrix.data(), eedf.data(), grid.cellNumber);
        } else {
            // TODO: Find a way to distinguish when to use LU and when to use Hessenberg reduction.

            // normal reduction wont work when vector contains first subdiagonal (if you want to use that
            // remove the first subdiagonal in evaluateMatrix().
//            LinAlg::hessenbergReductionOptimal(matrix.data(), &superElasticThresholds[0], grid.cellNumber,
//                                               superElasticThresholds.size());

//            Log<Message>::Notify("Is upper Hessenberg?");
//            Log<Message>::Notify(LinAlg::isUpperHessenberg(matrix.data(), grid.cellNumber));

//            eedf.setZero();
//            eedf[0] = 1.;

//            matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;

//            LinAlg::hessenberg(matrix.data(), eedf.data(), grid.cellNumber);

            auto *p = new uint32_t[grid.cellNumber];

            for (uint32_t i = 0; i < grid.cellNumber; ++i)
                p[i] = i;

            LinAlg::hessenbergReductionPartialPiv(matrix.data(), &superElasticThresholds[0], p, grid.cellNumber,
                                                  superElasticThresholds.size());

            eedf.setZero();
            eedf[0] = 1.;

            matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;

            LinAlg::hessenberg(matrix.data(), eedf.data(), p, grid.cellNumber);

            delete[] p;

//            Vector b = Vector::Zero(grid.cellNumber);
//            b[0] = 1;
//
//            matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;
//
//            eedf = matrix.partialPivLu().solve(b);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::cerr << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;

        eedf /= eedf.dot(grid.getCells().cwiseSqrt() * grid.step);
    }

    void ElectronKinetics::evaluateMatrix() {
        mixture.evaluateTotalAndElasticCS();

        evaluateElasticOperator();

        evaluateFieldOperator();

        if (!mixture.CARGasses.empty())
            evaluateCAROperator();

        evaluateInelasticOperators();

        evaluateIonizationOperator();

        evaluateAttachmentOperator();

        // Add the first subdiagonal to the vector.
//        superElasticThresholds.emplace_back(1);

        // Sort and erase duplicates.
        std::sort(superElasticThresholds.begin(), superElasticThresholds.end());
        superElasticThresholds.erase(unique(superElasticThresholds.begin(), superElasticThresholds.end()),
                                     superElasticThresholds.end());

        Log<Message>::Notify("Superelastic thresholds:");
        for (const auto entry : superElasticThresholds) {
            std::cerr << entry << std::endl;
        }
    }

    void ElectronKinetics::evaluateElasticOperator() {
        const double Tg = workingConditions->gasTemperature;

        const double factor1 = (Constant::kBeV * Tg / grid.step + 0.5) / grid.step;
        const double factor2 = (Constant::kBeV * Tg / grid.step - 0.5) / grid.step;

        g_c = grid.getNodes().cwiseAbs2().cwiseProduct(mixture.elasticCrossSection) * 2;

        g_c[0] = 0.;
        g_c[g_c.size() - 1] = 0.;

        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
            elasticMatrix(k, k) = -(g_c[k] * factor1 + g_c[k + 1] * factor2);

            if (k > 0)
                elasticMatrix(k, k - 1) = g_c[k] * factor2;

            if (k < grid.cellNumber - 1)
                elasticMatrix(k, k + 1) = g_c[k + 1] * factor1;
        }
    }

    void ElectronKinetics::evaluateFieldOperator() {
        const double EoN = workingConditions->reducedFieldSI;
        const double WoN = workingConditions->reducedExcFreqSI;
        const double me = Constant::electronMass;
        const double e = Constant::electronCharge;

        Vector &cs = mixture.totalCrossSection;

        g_E = ((EoN * EoN / 3) * grid.getNodes()).array() /
              (cs.array() + (me * WoN * WoN / (2 * e)) / (grid.getNodes().cwiseProduct(cs)).array());

        g_E[0] = 0.;
        g_E[g_E.size() - 1] = 0.;

        const double sqStep = grid.step * grid.step;

        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
            fieldMatrix(k, k) = -(g_E[k] + g_E[k + 1]) / sqStep;

            if (k > 0)
                fieldMatrix(k, k - 1) = g_E[k] / sqStep;

            if (k < grid.cellNumber - 1)
                fieldMatrix(k, k + 1) = g_E[k + 1] / sqStep;
        }
    }

    void ElectronKinetics::evaluateCAROperator() {
        const double Tg = workingConditions->gasTemperature;

        const double factor1 = (Constant::kBeV * Tg / grid.step + 0.5) / grid.step;
        const double factor2 = (Constant::kBeV * Tg / grid.step - 0.5) / grid.step;

        double sigma0B = 0.;

        for (auto *gas : mixture.CARGasses) {
            sigma0B += gas->fraction * gas->electricQuadrupoleMoment * gas->rotationalConstant;
        }

        sigma0B *= 8. * Constant::pi / (15. * Constant::electronCharge);
        g_CAR = grid.getNodes() * (4. * sigma0B);

        // Boundary conditions
        g_CAR[0] = 0.;
        g_CAR[grid.cellNumber] = 0.;

        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
            CARMatrix(k, k) = -(g_CAR[k] * factor1 + g_CAR[k + 1] * factor2);

            if (k > 0)
                CARMatrix(k, k - 1) = g_CAR[k] * factor2;

            if (k < grid.cellNumber - 1)
                CARMatrix(k, k + 1) = g_CAR[k + 1] * factor1;
        }
    }

    void ElectronKinetics::evaluateInelasticOperators() {
        const uint32_t cellNumber = grid.cellNumber;

        for (auto *gas : mixture.gasses) {
            for (auto vecIndex = (uint8_t) CollisionType::excitation;
                 vecIndex <= (uint8_t) CollisionType::rotational; ++vecIndex) {

                for (const auto *collision : gas->collisions[vecIndex]) {
                    const double threshold = collision->crossSection->threshold;

                    if (threshold < grid.step || threshold > grid.getNodes()[grid.cellNumber])
                        continue;

                    const double targetDensity = collision->target->density;

                    if (targetDensity != 0) {
                        const uint32_t numThreshold = std::floor(threshold / grid.step);

                        Vector cellCrossSection(cellNumber);

                        for (uint32_t i = 0; i < cellNumber; ++i)
                            cellCrossSection[i] =
                                    0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

                        for (uint32_t k = 0; k < cellNumber; ++k) {
                            if (k < cellNumber - numThreshold)
                                inelasticMatrix(k, k + numThreshold) +=
                                        targetDensity * grid.getCells()[k + numThreshold] *
                                        cellCrossSection[k + numThreshold];

                            inelasticMatrix(k, k) -= targetDensity * grid.getCells()[k] * cellCrossSection[k];
                        }

                        if (collision->isReverse) {
                            const double swRatio = collision->target->statisticalWeight /
                                                   collision->products[0]->statisticalWeight;
                            const double productDensity = collision->products[0]->density;

                            if (productDensity == 0)
                                continue;

                            if (numThreshold > 1)
                                superElasticThresholds.emplace_back(numThreshold);

                            if (numThreshold != 1) hasSuperelastics = true;

                            for (uint32_t k = 0; k < cellNumber; ++k) {
                                if (k >= numThreshold)
                                    inelasticMatrix(k, k - numThreshold) +=
                                            swRatio * productDensity * grid.getCells()[k] * cellCrossSection[k];

                                if (k < cellNumber - numThreshold)
                                    inelasticMatrix(k, k) -= swRatio * productDensity *
                                                             grid.getCells()[k + numThreshold] *
                                                             cellCrossSection[k + numThreshold];
                            }
                        }
                    }
                }
            }
        }
    }

    void ElectronKinetics::evaluateIonizationOperator() {
        bool hasValidCollisions = false;

        for (const auto *gas : mixture.gasses) {
            for (const auto *collision : gas->collisions[(uint8_t) CollisionType::ionization]) {
                const double threshold = collision->crossSection->threshold;

                if (threshold > grid.getNode(grid.cellNumber))
                    continue;

                hasValidCollisions = true;

                const double density = collision->target->density;
                const auto numThreshold = (uint32_t) (threshold / grid.step);

                Vector cellCrossSection(grid.cellNumber);

                for (uint32_t i = 0; i < grid.cellNumber; ++i)
                    cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

                switch (ionizationOperatorType) {
                    case IonizationOperatorType::conservative:
                        break;

                        // TODO: check this
                    case IonizationOperatorType::oneTakesAll:
                        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                            if (k < grid.cellNumber - numThreshold)
                                ionizationMatrix(k, k + numThreshold) +=
                                        density * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                            const double term = density * grid.getCell(k) * cellCrossSection(k);

                            ionizationMatrix(k, k) -= term;
                            ionizationMatrix(1, k) -= term;
                        }
                        break;

                        // TODO: check this
                    case IonizationOperatorType::equalSharing:
                        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                            ionizationMatrix(k, k) -= density * grid.getCell(k) * cellCrossSection[k];

                            if (k < (grid.cellNumber - numThreshold) / 2)
                                ionizationMatrix(k, 2 * k + numThreshold) +=
                                        4 * density * grid.getCell(2 * k + numThreshold) *
                                        cellCrossSection(2 * k + numThreshold);
                        }
                        break;
                    case IonizationOperatorType::sdcs:
                        double W = gas->OPBParameter;

                        if (W < 0)
                            W = threshold;

                        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                            const uint32_t end = std::min(2 * (k + 1) + numThreshold, grid.cellNumber);

                            if (k > numThreshold) {
                                const uint32_t half = (k + 1 - numThreshold) / 2;
                                const double numerator = 1 / std::atan((grid.getCell(k) - threshold) / (2 * W));
                                double sum = 0.;

                                for (uint32_t i = 0; i < half; ++i)
                                    sum += numerator / (W + grid.getCell(i) * grid.getCell(i) / W);

                                ionizationMatrix(k, k) -=
                                        density * grid.step * grid.getCell(k) * cellCrossSection[k] * sum;
                            }

                            if (k + numThreshold + 1 < grid.cellNumber) {
                                for (uint32_t i = k + numThreshold + 1; i < end; ++i) {
                                    ionizationMatrix(k, i) +=
                                            density * grid.step * grid.getCell(i) * cellCrossSection[i] /
                                            (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                             (W + std::pow(grid.getCell(i - k - numThreshold - 1), 2) / W));
                                }
                            }

                            // This last section might need some adjustments because of indexing
                            // differences between Matlab and C++ (since indexes are multiplied here).

                            for (uint32_t i = 2 * (k + 1) + numThreshold - 1; i < grid.cellNumber; ++i) {
                                ionizationMatrix(k, i) +=
                                        density * grid.step * grid.getCell(i) * cellCrossSection[i] /
                                        (std::atan((grid.getCell(i) - threshold) / (2 * W)) *
                                         (W + std::pow(grid.getCell(k), 2) / W));
                            }
                        }
                        break;
                }

                // Evaluation of the conservative ionization operator

                if (numThreshold == 0)
                    continue;

                for (uint32_t k = 0; k < grid.cellNumber; ++k) {
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

    void ElectronKinetics::evaluateAttachmentOperator() {
        attachmentMatrix.setZero();
        attachmentConservativeMatrix.setZero();

        const uint32_t cellNumber = grid.cellNumber;

        for (const auto *gas : mixture.gasses) {
            for (const auto *collision : gas->collisions[(uint8_t) CollisionType::attachment]) {
                const double threshold = collision->crossSection->threshold;

                if (threshold > grid.getNode(cellNumber))
                    continue;

                includeNonConservativeAttachment = true;

                const uint32_t numThreshold = std::floor(threshold / grid.step);

                Vector cellCrossSection(cellNumber);

                const double targetDensity = collision->target->density;

                for (uint32_t i = 0; i < cellNumber; ++i)
                    cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

                for (uint32_t k = 0; k < cellNumber; ++k)
                    attachmentMatrix(k, k) += targetDensity * grid.getCell(k) * cellCrossSection[k];

                if (numThreshold == 0) continue;

                for (uint32_t k = 0; k < cellNumber; ++k) {
                    if (k < cellNumber - numThreshold)
                        attachmentConservativeMatrix(k, k + numThreshold) =
                                targetDensity * grid.getCell(k + numThreshold) * cellCrossSection[k + numThreshold];

                    attachmentConservativeMatrix(k, k) = targetDensity * grid.getCell(k) * cellCrossSection[k];
                }
            }
        }
    }

    void ElectronKinetics::mixingDirectSolutions() {
        invertLinearMatrix();

        const uint32_t numCells = grid.cellNumber;

//        solveEEColl();
//        return;

        switch (growthModelType) {
            case GrowthModelType::spatial:
                ionSpatialGrowthD.setZero(numCells, numCells);
                ionSpatialGrowthU.setZero(numCells, numCells);
                fieldMatrixSpatGrowth.setZero(numCells, numCells);

                solveSpatialGrowthMatrix();
                break;

            case GrowthModelType::temporal:
                ionTemporalGrowth.setZero(numCells, numCells);
                fieldMatrixTempGrowth.setZero(numCells, numCells);

                solveTemporalGrowthMatrix();
                break;
        }
    }

    void ElectronKinetics::solveSpatialGrowthMatrix() {
        const double e = Constant::electronCharge,
                m = Constant::electronMass,
                EoN = workingConditions->reducedFieldSI;

        Vector cellTotalCrossSection(grid.cellNumber);

        for (uint32_t i = 0; i < grid.cellNumber; ++i)
            cellTotalCrossSection[i] = .5 * (mixture.totalCrossSection[i] + mixture.totalCrossSection[i + 1]);

        Matrix baseMatrix;

        if (!mixture.CARGasses.empty()) {
            baseMatrix = 1.e20 * (elasticMatrix + fieldMatrix + CARMatrix + inelasticMatrix + ionizationMatrix +
                                  attachmentMatrix);
        } else {
            baseMatrix = 1.e20 * (elasticMatrix + fieldMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
        }
        // TODO: add ee collision term here

        Vector baseDiag(grid.cellNumber), baseSubDiag(grid.cellNumber), baseSupDiag(grid.cellNumber);

        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
            baseDiag[k] = baseMatrix(k, k);

            if (k > 0)
                baseSubDiag[k] = baseMatrix(k, k - 1);

            if (k < grid.cellNumber - 1)
                baseSupDiag[k] = baseMatrix(k, k + 1);
        }

        Vector integrandCI = (sqrt(2. * e / m) * grid.step) * (ionizationMatrix + attachmentMatrix).colwise().sum();

        double CIEffNew = eedf.dot(integrandCI);
        double CIEffOld = CIEffNew / 3;

        CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

        // diffusion and mobility components of the spatial growth terms
        // can be removed since this is already done in the directMixing function
        ionSpatialGrowthD.setZero();
        ionSpatialGrowthU.setZero();

        Vector tempVector = grid.getCells().array() / (3. * cellTotalCrossSection).array();

        Vector D0 = tempVector,
                U0sup(grid.cellNumber),
                U0inf(grid.cellNumber);

        U0sup[0] = 0.;
        U0inf[grid.cellNumber - 1] = 0.;

        for (uint32_t j = 0; j < grid.cellNumber; ++j) {
            if (j != 0)
                U0sup[j] = EoN / (2. * grid.step) * tempVector[j - 1];

            if (j != grid.cellNumber - 1)
                U0inf[j] = -EoN / (2. * grid.step) * tempVector[j + 1];
        }

        Vector U0 = U0sup + U0inf;

        double ND = sqrt(2 * e / m) * grid.step * D0.dot(eedf),
                muE = -sqrt(2 * e / m) * grid.step * U0.dot(eedf);

        double alphaRedEffNew, alphaRedEffOld = 0.;

        if (muE * muE - 4 * CIEffNew * ND < 0.) {
            alphaRedEffNew = CIEffNew / muE;
        } else {
            alphaRedEffNew = (muE - sqrt(muE * muE - 4 * CIEffNew * ND)) / (2 * ND);
        }

        uint32_t iter = 0;
        bool hasConverged = false;

        const Vector g_fieldSpatialBase = (EoN / 6) * grid.getNodes().array() / mixture.totalCrossSection.array();

        while (!hasConverged) {
            g_fieldSpatialGrowth = alphaRedEffNew * g_fieldSpatialBase;
            g_fieldSpatialGrowth[0] = 0.;
            g_fieldSpatialGrowth[grid.cellNumber] = 0.;

            for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                fieldMatrixSpatGrowth(k, k) = (g_fieldSpatialGrowth[k + 1] - g_fieldSpatialGrowth[k]) / grid.step;

                if (k > 0)
                    fieldMatrixSpatGrowth(k, k - 1) = -g_fieldSpatialGrowth[k] / grid.step;

                if (k < grid.cellNumber - 1)
                    fieldMatrixSpatGrowth(k, k + 1) = g_fieldSpatialGrowth[k + 1] / grid.step;
            }

            for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                ionSpatialGrowthD(k, k) = alphaRedEffNew * alphaRedEffNew * D0[k];

                if (k > 0)
                    ionSpatialGrowthU(k, k - 1) = alphaRedEffNew * U0inf[k - 1];

                if (k < grid.cellNumber - 1)
                    ionSpatialGrowthU(k, k + 1) = alphaRedEffNew * U0sup[k + 1];
            }

            // TODO: add ee-col in following calculation

            for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                baseMatrix(k, k) = baseDiag[k] + 1.e20 * (fieldMatrixSpatGrowth(k, k) + ionSpatialGrowthD(k, k));

                if (k > 0)
                    baseMatrix(k, k - 1) =
                            baseSubDiag[k] + 1.e20 * (fieldMatrixSpatGrowth(k, k - 1) + ionSpatialGrowthU(k, k - 1));

                if (k < grid.cellNumber - 1)
                    baseMatrix(k, k + 1) =
                            baseSupDiag[k] + 1.e20 * (fieldMatrixSpatGrowth(k, k + 1) + ionSpatialGrowthU(k, k + 1));
            }

            Vector eedfNew = eedf;

            invertMatrix(baseMatrix);

            CIEffOld = CIEffNew;
            CIEffNew = eedf.dot(integrandCI);

            CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

            ND = sqrt(2 * e / m) * grid.step * D0.dot(eedf);
            muE = -sqrt(2 * e / m) * grid.step * U0.dot(eedf);

            const double discriminant = muE * muE - 4 * CIEffNew * ND;

            alphaRedEffOld = alphaRedEffNew;
            alphaRedEffNew = (discriminant < 0) ? CIEffNew / muE
                                                : (muE - sqrt(discriminant)) / (2 * ND);

            alphaRedEffNew = mixingParameter * alphaRedEffNew + (1 - mixingParameter) * alphaRedEffOld;

            if ((alphaRedEffNew == 0 || abs(alphaRedEffNew - alphaRedEffOld) / alphaRedEffOld < 1.e-10) &&
                (((eedf - eedfNew).cwiseAbs().array() / eedf.array()).maxCoeff() < maxEedfRelError || iter > 150)) {
                hasConverged = true;

                if (iter > 150 && !includeEECollisions)
                    Log<Message>::Warning("Iterative spatial growth scheme did not converge.");
            }

            ++iter;
        }

        std::cerr << "Number of iterations to convergence: " << iter << ".\n";

        alphaRedEff = alphaRedEffOld;
        CIEff = CIEffOld;
    }

    void ElectronKinetics::solveTemporalGrowthMatrix() {
        const double e = Constant::electronCharge,
                m = Constant::electronMass,
                EoN = workingConditions->reducedFieldSI,
                WoN = workingConditions->reducedExcFreqSI;

        Matrix baseMatrix;

        if (!mixture.CARGasses.empty()) {
            baseMatrix = 1.e20 * (elasticMatrix + CARMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
        } else {
            baseMatrix = 1.e20 * (elasticMatrix + inelasticMatrix + ionizationMatrix + attachmentMatrix);
        }

        Vector baseDiag(grid.cellNumber), baseSubDiag(grid.cellNumber), baseSupDiag(grid.cellNumber);

        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
            baseDiag[k] = baseMatrix(k, k);

            if (k > 0)
                baseSubDiag[k] = baseMatrix(k, k - 1);

            if (k < grid.cellNumber - 1)
                baseSupDiag[k] = baseMatrix(k, k + 1);
        }

        Vector integrandCI = (sqrt(2. * e / m) * grid.step) * (ionizationMatrix + attachmentMatrix).colwise().sum();

        double CIEffNew = eedf.dot(integrandCI);
        double CIEffOld = CIEffNew / 3.;

        CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

        Vector totalCSI(grid.cellNumber + 1),
                eedfNew(grid.cellNumber);

        bool hasConverged = false;
        uint32_t iter = 0;

        while (!hasConverged) {

            totalCSI[0] = mixture.totalCrossSection[0];

            const long double growthFactor = CIEffNew * sqrt(m / (2 * e));

            for (uint32_t i = 1; i <= grid.cellNumber; ++i) {
                totalCSI[i] = mixture.totalCrossSection[i] + growthFactor / sqrt(i * grid.step);
            }

            g_fieldTemporalGrowth = ((EoN * EoN / 3) * grid.getNodes()).array() / (totalCSI.array() +
                                                                                   (m * WoN * WoN / (2 * e)) /
                                                                                   (grid.getNodes().cwiseProduct(
                                                                                           totalCSI)).array());
            g_fieldTemporalGrowth[0] = 0.;
            g_fieldTemporalGrowth[grid.cellNumber] = 0.;

            const double sqrStep = grid.step * grid.step;

            for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                fieldMatrixTempGrowth(k, k) = -(g_fieldTemporalGrowth[k] + g_fieldTemporalGrowth[k + 1]) / sqrStep;

                if (k > 0)
                    fieldMatrixTempGrowth(k, k - 1) = g_fieldTemporalGrowth[k] / sqrStep;

                if (k < grid.cellNumber - 1)
                    fieldMatrixTempGrowth(k, k + 1) = g_fieldTemporalGrowth[k + 1] / sqrStep;

                ionTemporalGrowth(k, k) = -growthFactor * sqrt(grid.getCell(k));
            }

            for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                baseMatrix(k, k) = baseDiag[k] + 1.e20 * (fieldMatrixTempGrowth(k, k) + ionTemporalGrowth(k, k));

                if (k > 0)
                    baseMatrix(k, k - 1) = baseSubDiag[k] + 1.e20 * fieldMatrixTempGrowth(k, k - 1);

                if (k < grid.cellNumber - 1)
                    baseMatrix(k, k + 1) = baseSupDiag[k] + 1.e20 * fieldMatrixTempGrowth(k, k + 1);
            }

            eedfNew = eedf;

            invertMatrix(baseMatrix);

            CIEffOld = CIEffNew;
            CIEffNew = eedf.dot(integrandCI);
            CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

            if ((CIEffNew == 0 || abs(CIEffNew - CIEffOld) / CIEffOld < 1.e10) &&
                (((eedf - eedfNew).cwiseAbs().array() / eedf.array()).maxCoeff() < maxEedfRelError || iter > 150)) {
                hasConverged = true;

                if (iter > 150 && !includeEECollisions)
                    Log<Message>::Warning("Iterative temporal growth scheme did not converge.");
            }

            ++iter;
        }

        std::cerr << "Number of iterations to convergence: " << iter << ".\n";

        CIEff = CIEffOld;
    }

    void ElectronKinetics::solveEEColl() {
        const double e = Constant::electronCharge,
                e0 = Constant::vacuumPermittivity,
                ne = workingConditions->electronDensity,
                n0 = workingConditions->gasDensity,
                EoN = workingConditions->reducedField;

        Matrix Mee = Matrix::Zero(grid.cellNumber, grid.cellNumber);

        Matrix baseMatrix;

        // Splitting all possible options for the best performance.

        if (includeNonConservativeIonization || includeNonConservativeAttachment) {
            if (growthModelType == GrowthModelType::spatial) {
                if (mixture.CARGasses.empty()) {
                    baseMatrix = 1.e20 * (ionizationMatrix + attachmentMatrix + elasticMatrix + inelasticMatrix +
                                          fieldMatrix + ionSpatialGrowthD + ionSpatialGrowthU + fieldMatrixSpatGrowth);
                } else {
                    baseMatrix =
                            1.e20 * (ionizationMatrix + attachmentMatrix + elasticMatrix + inelasticMatrix + CARMatrix +
                                     fieldMatrix + ionSpatialGrowthD + ionSpatialGrowthU + fieldMatrixSpatGrowth);
                }
            } else if (growthModelType == GrowthModelType::temporal) {
                if (mixture.CARGasses.empty()) {
                    baseMatrix = 1.e20 * (ionizationMatrix + attachmentMatrix + elasticMatrix + inelasticMatrix +
                                          fieldMatrix + ionTemporalGrowth + fieldMatrixTempGrowth);
                } else {
                    baseMatrix =
                            1.e20 * (ionizationMatrix + attachmentMatrix + elasticMatrix + inelasticMatrix + CARMatrix +
                                     fieldMatrix + ionTemporalGrowth + fieldMatrixTempGrowth);
                }
            }
        } else {
            if (mixture.CARGasses.empty()) {
                baseMatrix = 1.e20 * (ionConservativeMatrix + attachmentConservativeMatrix + elasticMatrix +
                                      inelasticMatrix + fieldMatrix);
            } else {
                baseMatrix = 1.e20 * (ionConservativeMatrix + attachmentConservativeMatrix + elasticMatrix +
                                      inelasticMatrix + fieldMatrix + CARMatrix);
            }
        }

        // Storing the initial diagonals of the matrix in three separate vectors.
        // This allows us to skip the usage of 'matrixAux', saving a good amount of
        // memory for bigger matrices.

        Vector baseDiag(grid.cellNumber), baseSubDiag(grid.cellNumber), baseSupDiag(grid.cellNumber);

        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
            baseDiag[k] = baseMatrix(k, k);

            if (k > 0)
                baseSubDiag[k] = baseMatrix(k, k - 1);

            if (k < grid.cellNumber - 1)
                baseSupDiag[k] = baseMatrix(k, k + 1);
        }

        Bee.setZero(grid.cellNumber, grid.cellNumber);

        const Vector cellsThreeOverTwo = grid.getCells().cwiseProduct(grid.getCells().cwiseSqrt()),
                energyArray = -(grid.step / 2.) * grid.getCells().cwiseSqrt() + (2. / 3.) * cellsThreeOverTwo;

        for (uint32_t j = 0; j < grid.cellNumber - 1; ++j) {

            for (uint32_t i = 1; i <= j; ++i)
                Bee(i, j) = energyArray[i];

            const double value = 2. / 3. * std::pow(grid.getNode(j + 1), 1.5);

            for (uint32_t i = j + 1; i < grid.cellNumber; ++i)
                Bee(i, j) = value;
        }

        // detailed balance condition

        for (uint32_t j = 0; j < grid.cellNumber - 1; ++j) {
            for (uint32_t i = 1; i < grid.cellNumber; ++i) {
                Bee(i, j) = sqrt(Bee(i, j) * Bee(j + 1, i - 1));
            }
        }

        Aee = Bee.transpose();

        double meanEnergy = grid.step * cellsThreeOverTwo.dot(eedf),
                Te = 2. / 3. * meanEnergy,
                logC = std::log(12 * Constant::pi * std::pow(e0 * Te / e, 1.5) / std::sqrt(ne)),
                alpha = (ne / n0) * (e * e / (8 * Constant::pi * e0 * e0)) * logC;

        double ratioNew = 0.;
        Vector eedfNew = eedf;

        bool hasConverged = false;
        uint32_t iter = 0;

        Vector MeeDiag(grid.cellNumber), MeeSub(grid.cellNumber), MeeSup(grid.cellNumber);

        // In this implementation we completely skip the Mee matrix, saving both memory and time.

        while (!hasConverged) {
            Vector A = (alpha / grid.step) * (Aee * eedf),
                    B = (alpha / grid.step) * (Bee * eedf);

            for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                baseMatrix(k, k) = baseDiag[k] - 1.e20 * (A[k] + B[k]);

                if (k > 0)
                    baseMatrix(k, k - 1) = baseSubDiag[k] + 1.e20 * A[k - 1];

                if (k < grid.cellNumber - 1)
                    baseMatrix(k, k + 1) = baseSupDiag[k] + 1.e20 * B[k + 1];
            }

            invertMatrix(baseMatrix);

            evaluatePower(false);

            const double ratio = std::abs(power.electronElectron / power.reference);

            if (((eedf - eedfNew).cwiseAbs().array() / eedf.array()).maxCoeff() < maxEedfRelError) {
                if (std::abs(ratio) < 1.e-9) {
                    hasConverged = true;
                } else if (std::abs(ratio) > 1.e-9 && iter > 200) {
                    // TODO: give proper warnings
                    std::cerr << "error in ee-col convergence" << std::endl;
                    hasConverged = true;
                }
            } else if (iter == 300 && !(includeNonConservativeAttachment || includeNonConservativeIonization)) {
                hasConverged = true;
                Log<Message>::Warning("Electron-electron iterative scheme did not converge.");
            }

            if (iter > 0 && hasConverged) {
                double ratioOld = ratioNew;
                ratioNew = ratio;

                Vector eedfOld = eedfNew;
                eedfNew = eedf;

                const double norm = eedf.dot(grid.getCells().cwiseSqrt()) * grid.step;

                eedf = eedfNew - (ratioNew / (ratioNew - ratioOld)) * (eedfNew - eedfOld);

                for (uint32_t i = 0; i < grid.cellNumber; ++i) {
                    if (std::abs(norm * eedf[i]) < std::numeric_limits<double>::epsilon())
                        eedf[i] = std::abs(eedf[i]);
                        // TODO: ask if this is correct (might be an error in Matlab)
                    else if (eedf[i] < 0)
                        eedf[i] = abs(eedf[i]);
                }
            }

            meanEnergy = grid.step * cellsThreeOverTwo.dot(eedf);
            Te = 2. / 3. * meanEnergy;
            logC = std::log(12 * Constant::pi * std::pow(e0 * Te / e, 1.5) / std::sqrt(ne));
            alpha = (ne / n0) * (e * e / (8 * Constant::pi * e0 * e0)) * logC;

            iter++;
        }

        // TODO: save alpha
    }

    void ElectronKinetics::evaluatePower(bool isFinalSolution) {
        const double factor = sqrt(2. * Constant::electronCharge / Constant::electronMass),
                kTg = Constant::kBeV * workingConditions->gasTemperature,
                auxHigh = kTg + grid.step * .5,
                auxLow = kTg + grid.step * .5;

        double elasticNet = 0., elasticGain = 0.;

        for (uint32_t k = 0; k < grid.cellNumber; ++k) {
            elasticNet += eedf[k] * (g_c[k + 1] * auxLow - g_c[k] * auxLow);
            elasticGain += eedf[k] * (g_c[k + 1] - g_c[k]);
        }

        power.elasticNet = factor * elasticNet;
        power.elasticGain = factor * kTg * elasticGain;
        power.elasticLoss = power.elasticNet - power.elasticGain;

        if (!mixture.CARGasses.empty()) {
            double carNet = 0., carGain = 0.;

            for (uint32_t k = 0; k < grid.cellNumber - 1; ++k) {
                carNet += eedf[k] * (g_CAR[k + 1] * auxLow - g_CAR[k] * auxHigh);
                carGain += eedf[k] * (g_CAR[k + 1] - g_CAR[k]);
            }

            power.carNet = factor * carNet;
            power.carGain = factor * kTg * carGain;
            power.carLoss = power.carNet - power.carGain;
        }

        if (includeNonConservativeIonization || includeNonConservativeAttachment) {
            if (growthModelType == GrowthModelType::temporal) {
                double field = 0., growthModel = 0.;

                for (uint32_t k = 0; k < grid.cellNumber - 1; ++k) {
                    field += eedf[k] * (g_fieldTemporalGrowth[k + 1] - g_fieldTemporalGrowth[k]);
                    growthModel += eedf[k] * grid.getCell(k) * sqrt(grid.getCell(k));
                }

                power.field = factor * field;
                power.eDensGrowth = CIEff * grid.step * growthModel;
            } else if (growthModelType == GrowthModelType::spatial) {
                double field = 0., correction = 0., powerDiffusion = 0., powerMobility = 0.;
                Vector cellCrossSection(grid.cellNumber);

                for (uint32_t k = 0; k < grid.cellNumber - 1; ++k) {
                    field += eedf[k] * (g_E[k + 1] - g_E[k]);
                    correction -= eedf[k] * (g_fieldSpatialGrowth[k + 1] * auxLow -
                                             g_fieldSpatialGrowth[k] * auxLow);

                    // Diffusion and Mobility contributions
                    cellCrossSection[k] = .5 * (mixture.totalCrossSection[k] + mixture.totalCrossSection[k + 1]);
                    powerDiffusion += grid.getCell(k) * grid.getCell(k) * eedf[k] / cellCrossSection[k];

                    if (k > 0) {
                        powerMobility +=
                                grid.getCell(k) * grid.getCell(k) * (eedf[k + 1] - eedf[k - 1]) / cellCrossSection[k];
                    }
                }

                power.field = factor * (field + grid.step * correction);
                power.eDensGrowth = alphaRedEff * alphaRedEff * factor * grid.step / 3. * powerDiffusion +
                                    factor * alphaRedEff * (workingConditions->reducedFieldSI / 6.) *
                                    (grid.getCell(0) * grid.getCell(0) * eedf[1] / cellCrossSection[0] -
                                     grid.getCell(grid.cellNumber - 1) * grid.getCell(grid.cellNumber - 1) *
                                     eedf[grid.cellNumber - 2] / cellCrossSection[grid.cellNumber - 1] + powerMobility);
            }
        } else {
            double field = 0;

            for (uint32_t k = 0; k < grid.cellNumber - 1; ++k) {
                field += eedf[k] * (g_E[k + 1] - g_E[k]);
            }

            power.field = factor * field;
        }

        if (includeEECollisions) {
            power.electronElectron = (-factor * grid.step * grid.step) * ((Aee - Bee) * eedf).sum();
        }

        // Evaluate power absorbed per electron at unit gas density due to in- and superelastic collisions.
        for (auto *gas : mixture.gasses) {
            gas->evaluatePower(ionizationOperatorType, eedf);
            power += gas->getPower();
        }

        // TODO: The power per gas is not stored in the main power structure for now.

        auto *powerPtr = (double *) &power;

        double totalGain = 0., totalLoss = 0.;

        // Loop over the first 23 double member variables of the power struct.
        for (uint32_t i = 0; i < 22; ++i) {
            if (powerPtr[i] > 0)
                totalGain += powerPtr[i];
            else
                totalLoss += powerPtr[i];
        }

        power.balance = power.field + power.elasticNet + power.carNet + power.inelastic + power.superelastic +
                        power.eDensGrowth + power.electronElectron;
        power.relativeBalance = abs(power.balance) / totalGain;
        power.reference = totalGain;

        if (isFinalSolution && power.relativeBalance > maxPowerBalanceRelError)
            Log<PowerBalanceError>::Warning(maxPowerBalanceRelError);
    }

    void ElectronKinetics::plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
                                const Vector &x, const Vector &y) {
        std::cout << "unset key" << std::endl;
        std::cout << "set xlabel \"" << xlabel << "\"" << std::endl;
        std::cout << "set ylabel \"" << ylabel << "\"" << std::endl;
        std::cout << "set title \"" << title << "\"" << std::endl;
        std::cout << "set xrange [" << x[0] << ":" << x[x.size() - 1] << "]" << std::endl;
        // std::cout << "set logscale y" << std::endl;
        std::cout << "plot '-' w l" << std::endl;
        for (uint32_t i = 0; i < x.size(); ++i) {
            std::cout << x[i] << "\t" << y[i] << '\n';
        }
        std::cout << "e" << std::endl;
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