//
// Created by daan on 13-5-19.
//

#include "ElectronKinetics.h"
#include <chrono>
#include <cmath>

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
        matrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;

        auto begin = std::chrono::high_resolution_clock::now();

        if (!hasSuperelastics && !includeEECollisions) {
            eedf.setZero();
            eedf[0] = 1.;

            LinAlg::hessenberg(matrix.data(), eedf.data(), grid.cellNumber);
        } else {
            Vector b = Vector::Zero(grid.cellNumber);
            b[0] = 1;

            eedf = matrix.partialPivLu().solve(b);
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

                            hasSuperelastics = true;

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

        switch (growthModelType) {
            case GrowthModelType::spatial:
                ionSpatialGrowthD.setZero(numCells, numCells);
                ionSpatialGrowthU.setZero(numCells, numCells);
                fieldMatrixSpatGrowth.setZero(numCells, numCells);

                solveSpatialGrowthMatrix();
                break;

            case GrowthModelType::temporal:
                Log<Message>::Error("Temporal growth is not yet supported.");
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

        Vector integrandCI = (sqrt(2. * e / m) * grid.step) * (ionizationMatrix).colwise().sum();

        double CIEffNew = eedf.dot(integrandCI);
        double CIEffOld = CIEffNew / 3;

        CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

        // diffusion and mobility components of the spatial growth terms
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

        double alphaRedEffNew;

        if (muE * muE - 4 * CIEffNew * ND < 0.) {
            alphaRedEffNew = CIEffNew / muE;
        } else {
            alphaRedEffNew = (muE - sqrt(muE * muE - 4 * CIEffNew * ND)) / (2 * ND);
        }

//        printf("CIEffNew = %.16e\n", CIEffNew);
//        printf("alphaRedEffNew = %.16e\n", alphaRedEffNew);

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

//            Matrix boltzmannMatrix = baseMatrix;

            // TODO: add ee-col in following calculation

            for (uint32_t k = 0; k < grid.cellNumber; ++k) {
                baseMatrix(k, k) = baseDiag[k] + 1.e20 * (fieldMatrixSpatGrowth(k, k) + ionSpatialGrowthD(k, k));
//                boltzmannMatrix(k, k) += 1.e20 * (fieldMatrixSpatGrowth(k, k) + ionSpatialGrowthD(k, k));

                if (k > 0)
                    baseMatrix(k, k - 1) =
                            baseSubDiag[k] + 1.e20 * (fieldMatrixSpatGrowth(k, k - 1) + ionSpatialGrowthU(k, k - 1));
//                    boltzmannMatrix(k, k - 1) +=
//                            1.e20 * (fieldMatrixSpatGrowth(k, k - 1) + ionSpatialGrowthU(k, k - 1));

                if (k < grid.cellNumber - 1)
                    baseMatrix(k, k + 1) =
                            baseSupDiag[k] + 1.e20 * (fieldMatrixSpatGrowth(k, k + 1) + ionSpatialGrowthU(k, k +  1));
//                    boltzmannMatrix(k, k + 1) +=
//                            1.e20 * (fieldMatrixSpatGrowth(k, k + 1) + ionSpatialGrowthU(k, k + 1));
            }

            Vector eedfNew = eedf;

//            invertMatrix(boltzmannMatrix);
            invertMatrix(baseMatrix);

//            if (iter == 16) {
//                for (uint32_t i = 0; i < grid.cellNumber; ++i) {
//                    printf("%.16e\n", eedf[i]);
//                }
//            }

//            if (iter == 16) {
//                for (uint32_t i = 0; i < grid.cellNumber; ++i) {
//                    for (uint32_t j = 0; j < grid.cellNumber; ++j) {
//                        printf("%.16e", boltzmannMatrix(i, j));
//
//                        if (j != grid.cellNumber - 1)
//                            std::cout << '\t';
//                    }
//
//                    std::cout << '\n';
//                }
//            }
//            return;

            CIEffOld = CIEffNew;
            CIEffNew = eedf.dot(integrandCI);

            CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;

            ND = sqrt(2 * e / m) * grid.step * D0.dot(eedf);
            muE = -sqrt(2 * e / m) * grid.step * U0.dot(eedf);

            const double discriminant = muE * muE - 4 * CIEffNew * ND;

            const double alphaRedEffOld = alphaRedEffNew;
            alphaRedEffNew = (discriminant < 0) ? CIEffNew / muE
                                                : (muE - sqrt(discriminant)) / (2 * ND);

            alphaRedEffNew = mixingParameter * alphaRedEffNew + (1 - mixingParameter) * alphaRedEffOld;

//            printf("CIEffNew = %.16e\n", CIEffNew);
//            printf("alphaRedEffNew = %.16e\n", alphaRedEffNew);

            if ((alphaRedEffNew == 0 || abs(alphaRedEffNew - alphaRedEffOld) / alphaRedEffOld < 1.e-10) &&
                (((eedf - eedfNew).cwiseAbs().array() / eedf.array()).maxCoeff() < maxEedfRelError || iter > 150)) {
                hasConverged = true;

                if (iter > 150 && !includeEECollisions)
                    Log<Message>::Warning("Iterative spatial growth scheme did not converge.");
            }

            ++iter;
        }

        std::cerr << "Number of iterations to convergence: " << iter << ".\n";
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

/*
Vector totalCSI(grid.cellNumber + 1);

    while (!hasConverged)
    {

        totalCSI[0] = mixture.totalCrossSection[0];

        const double factor = CIEffNew * sqrt(m / (2 * e * grid.step));

        for (uint32_t i = 1; i <= grid.cellNumber; ++i)
        {
            totalCSI[i] = mixture.totalCrossSection[i] * (factor / sqrt(i));
        }
    } 
    */