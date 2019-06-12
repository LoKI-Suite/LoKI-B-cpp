//
// Created by daan on 13-5-19.
//

#include "ElectronKinetics.h"
#include <chrono>
#include <iomanip>

namespace loki {
    ElectronKinetics::ElectronKinetics(const ElectronKineticsSetup &setup, const WorkingConditions *workingConditions)
            : workingConditions(workingConditions), grid(setup.numerics.energyGrid), mixture(&grid),
              g_c(grid.cellNumber), eedf(grid.cellNumber), elasticMatrix(grid.cellNumber, grid.cellNumber),
              continuousMatrix(grid.cellNumber, grid.cellNumber) {

        mixture.initialize(setup, workingConditions);

        this->eedfType = setup.eedfType;
        this->shapeParameter = setup.shapeParameter;
        this->ionizationOperatorType = setup.ionizationOperatorType;
        this->growthModelType = setup.growthModelType;
        this->includeEECollisions = setup.includeEECollisions;

//        this->plot("Total Elastic Cross Section N2", "Energy (eV)", "Cross Section (m^2)",
//                   mixture.grid->getNodes(), mixture.totalCrossSection);

        elasticMatrix.setZero(grid.cellNumber, grid.cellNumber);
        fieldMatrix.setZero(grid.cellNumber, grid.cellNumber);
        inelasticMatrix.setZero(grid.cellNumber, grid.cellNumber);

        if (!mixture.CARGasses.empty())
            CARMatrix.setZero(grid.cellNumber, grid.cellNumber);

        this->evaluateMatrix();
    }

    void ElectronKinetics::solve() {
        Matrix boltzmannMatrix = (elasticMatrix + fieldMatrix + CARMatrix + inelasticMatrix).array() * 1.e20;
        Vector b = Vector::Zero(grid.cellNumber);

        // Induce normalization condition

        boltzmannMatrix.row(0) = grid.getCells().cwiseSqrt() * grid.step;
        b[0] = 1;

        auto begin = std::chrono::high_resolution_clock::now();

        eedf = boltzmannMatrix.partialPivLu().solve(b);

        auto end = std::chrono::high_resolution_clock::now();
        std::cerr << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;

        eedf /= eedf.dot(grid.getCells().cwiseSqrt() * grid.step);

//        for (uint32_t i = 0; i < eedf.size(); ++i) {
//            printf("%.16e\n", eedf[i]);
//        }

        this->plot("Eedf due to elastic collisions", "Energy (eV)", "Eedf (Au)", grid.getCells(), eedf);
    }

    void ElectronKinetics::evaluateMatrix() {
        mixture.evaluateTotalAndElasticCS();

        evaluateElasticOperator();

        evaluateFieldOperator();

        if (!mixture.CARGasses.empty())
            evaluateCAROperator();

        evaluateInelasticOperators();
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

                    const uint32_t numThreshold = std::floor(threshold / grid.step);

                    Vector cellCrossSection(cellNumber);

                    const double targetDensity = collision->target->density;

                    for (uint32_t i = 0; i < cellNumber; ++i)
                        cellCrossSection[i] = 0.5 * ((*collision->crossSection)[i] + (*collision->crossSection)[i + 1]);

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

    void ElectronKinetics::plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
                                const Vector &x, const Vector &y) {
        std::cout << "unset key" << std::endl;
        std::cout << "set xlabel \"" << xlabel << "\"" << std::endl;
        std::cout << "set ylabel \"" << ylabel << "\"" << std::endl;
        std::cout << "set title \"" << title << "\"" << std::endl;
        std::cout << "set xrange [" << x[0] << ":" << x[x.size() - 1] << "]" << std::endl;
//        std::cout << "set logscale y" << std::endl;
        std::cout << "plot '-' w l" << std::endl;
        for (uint32_t i = 0; i < x.size(); ++i) {
            std::cout << x[i] << "\t" << y[i] << '\n';
        }
        std::cout << "e" << std::endl;
    }
}