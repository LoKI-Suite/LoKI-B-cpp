//
// Created by daan on 21-5-19.
//

#include "EedfCollision.h"
#include "Constant.h"
#include "Log.h"

namespace loki {

    EedfCollision::EedfCollision(Enumeration::CollisionType type, std::vector<EedfState *> &reactants,
                                 std::vector<EedfState *> &products, std::vector<uint16_t> &stoiCoeff, bool isReverse)
            : Collision(type, reactants[0], products, stoiCoeff, isReverse) {

        if (reactants.size() != 1)
            Log<MultipleReactantInEedfCol>::Warning(*this);

        if (isReverse && products.size() > 1)
            Log<MultipleProductsInReverse>::Error(*this);
    }

    EedfCollision::~EedfCollision() {
        delete crossSection;
    }

    bool EedfCollision::operator==(const EedfCollision &other) {
        if (type != other.type) return false;
        if (target != other.target) return false;

        // Comparing pointers since there is only one instantiation of each state.
        for (const auto *state : other.products) {
            if (std::find(products.begin(), products.end(), state) == products.end())
                return false;
        }

        return true;
    }

    EedfState *EedfCollision::getTarget() {
        return target;
    }

    std::ostream &operator<<(std::ostream &os, const EedfCollision &collision) {
        os << "e + " << *collision.target
           << (collision.isReverse ? " <->" : " ->");

        if (collision.type != CollisionType::attachment) os << " e +";
        if (collision.type == CollisionType::ionization) os << " e +";

        for (uint32_t i = 0; i < collision.products.size(); ++i) {
            os << ' ' << *collision.products[i];
            if (i < collision.products.size() - 1) os << " +";
        }

        return os;
    }

    void EedfCollision::superElastic(const Vector &energyData, Vector &result) const {
        if (!isReverse) {
            Log<SuperElasticForNonReverse>::Error(*this);
        }

        result.resize(energyData.size());

        Vector superElasticEnergies = energyData.array() + crossSection->threshold;

        crossSection->interpolate(superElasticEnergies, result);

        const double swRatio = target->statisticalWeight / products[0]->statisticalWeight;

        uint32_t minIndex = 0;

        if (energyData[0] == 0)
            ++minIndex;

        for (uint32_t i = minIndex; i < result.size(); ++i) {
            result[i] *= swRatio * (1 + (crossSection->threshold / energyData[i]));
        }

        if (energyData[0] == 0) result[0] = 0;
    }

    CollPower EedfCollision::evaluateConservativePower(const Vector &eedf) {
        const Grid *grid = crossSection->getGrid();
        const uint32_t n = grid->cellNumber;

        CollPower collPower;

        Vector cellCrossSection(n);

        for (uint32_t i = 0; i < n; ++i) {
            cellCrossSection[i] = .5 * ((*crossSection)[i] + (*crossSection)[i + 1]);
        }

        auto lmin = (uint32_t) (crossSection->threshold / grid->step);

        const double factor = sqrt(2 * Constant::electronCharge / Constant::electronMass);

        double ineSum = 0;

        for (uint32_t i = lmin; i < n; ++i) {
            ineSum += eedf[i] * grid->getCell(i) * cellCrossSection[i];
        }

        collPower.ine = -factor * target->density * grid->step * grid->getNode(lmin) * ineSum;

        if (isReverse) {
            const double statWeightRatio = target->statisticalWeight / products[0]->statisticalWeight;

            double supSum = 0;

            for (uint32_t i = lmin; i < n; ++i) {
                supSum += eedf[i - lmin] * grid->getCell(i) * cellCrossSection[i];
            }

            collPower.sup +=
                    factor * statWeightRatio * products[0]->density * grid->step * grid->getNode(lmin) * supSum;
        }

        return collPower;
    }

    CollPower EedfCollision::evaluateNonConservativePower(const Vector &eedf,
                                                          const IonizationOperatorType ionizationOperatorType,
                                                          const double OPBParameter) {
        const Grid *grid = crossSection->getGrid();
        const uint32_t n = grid->cellNumber;

        CollPower collPower;

        const double factor = sqrt(2 * Constant::electronCharge / Constant::electronMass);

        Vector cellCrossSection(n);

        for (uint32_t i = 0; i < n; ++i) {
            cellCrossSection[i] = .5 * ((*crossSection)[i] + (*crossSection)[i + 1]);
        }

        auto lmin = (uint32_t) (crossSection->threshold / grid->step);

        if (type == CollisionType::ionization) {

            if (ionizationOperatorType == IonizationOperatorType::equalSharing) {
                double sumOne = 0., sumTwo = 0., sumThree = 0.;

                for (uint32_t i = lmin - 1; i < n; ++i) {
                    sumOne += grid->getCell(i) * grid->getCell(i) * cellCrossSection[i] * eedf[i];
                }

                for (uint32_t i = 1 + lmin; i < n; i += 2) {
                    const double term = grid->getCell(i) * cellCrossSection[i] * eedf[i];

                    sumTwo += term;
                    sumThree += grid->getCell(i) * term;
                }

                collPower.ine = -factor * target->density * grid->step *
                                (sumOne + 2 * grid->getCell(lmin) * sumTwo - 2 * sumThree);

            } else if (ionizationOperatorType == IonizationOperatorType::oneTakesAll) {
                double sum = 0.;

                for (uint32_t i = lmin - 1; i < n; ++i) {
                    sum += grid->getCell(i) * cellCrossSection[i] * eedf[i];
                }

                collPower.ine = -factor * target->density * grid->step * grid->getCell(lmin - 1) * sum;
            } else if (ionizationOperatorType == IonizationOperatorType::sdcs) {
                double w = OPBParameter;

                if (w == 0) w = crossSection->threshold;

                Vector TICS = Vector::Zero(grid->cellNumber);
                Vector auxVec = 1. / (1. + (grid->getCells().cwiseAbs2().array() / (w * w)));

                for (int64_t k = 1; k < n; ++k) {
                    auxVec[k] = auxVec[k] + auxVec[k - 1];
                    int64_t kmax = (k + 1 - lmin) / 2;

                    if (kmax > 0)
                        TICS[k] += cellCrossSection[k] * auxVec[kmax-1] /
                                   (w * atan((grid->getCell(k) - crossSection->threshold) / (2 * w)));
                }

                collPower.ine = -factor * target->density * grid->getCell(lmin) * grid->step *
                                eedf.cwiseProduct(grid->getCells().cwiseProduct(grid->step * TICS)).sum();
            }
        } else if (type == CollisionType::attachment) {
            double sum = 0.;

            for (uint32_t i = lmin; i < n; ++i) {
                sum += eedf[i] * grid->getCell(i) * grid->getCell(i) * cellCrossSection[i];
            }

            collPower.ine = -factor * target->density * grid->step * sum;
        }// else {
//             error
//        }

        return collPower;
    }

    RateCoefficient EedfCollision::evaluateRateCoefficient(const Vector &eedf) {
        const double factor = std::sqrt(2. * Constant::electronCharge / Constant::electronMass);
        const Grid *grid = crossSection->getGrid();

        const uint32_t nNodes = grid->cellNumber + 1,
                nCells = grid->cellNumber;

        const auto lmin = (uint32_t) (crossSection->threshold / grid->step);

        Vector cellCrossSection =
                .5 * (crossSection->segment(lmin, nNodes - 1 - lmin) + crossSection->tail(nNodes - 1 - lmin));

        ineRateCoeff = factor * grid->step * cellCrossSection.cwiseProduct(grid->getCells().tail(nCells - lmin)).dot(
                eedf.tail(nCells - lmin));

        if (isReverse) {
            const double tStatWeight = target->statisticalWeight,
                    pStatWeight = products[0]->statisticalWeight;

            //TODO: check if one of the two is zero

            const double statWeightRatio = tStatWeight / pStatWeight;

            supRateCoeff = factor * statWeightRatio * grid->step *
                          cellCrossSection.cwiseProduct(grid->getCells().tail(nCells - lmin)).dot(
                                  eedf.head(nCells - lmin));
        }

        return {this, ineRateCoeff, supRateCoeff};
    }
}