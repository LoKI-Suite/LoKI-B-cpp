//
// Created by daan on 21-5-19.
//

#include "LoKI-B/EedfCollision.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Log.h"

namespace loki {

    EedfCollision::EedfCollision(CollisionType type,
        const StateVector &lhsStates,
        const CoeffVector &lhsCoeffs,
        const StateVector &rhsStates,
        const CoeffVector &rhsCoeffs,
        bool isReverse)
    : Collision(type, lhsStates, lhsCoeffs, rhsStates, rhsCoeffs, isReverse),
    m_lhsHeavyStates(lhsStates),
    m_lhsHeavyCoeffs(lhsCoeffs),
    m_rhsHeavyStates(rhsStates),
    m_rhsHeavyCoeffs(rhsCoeffs)
    {
        if (lhsStates.size() != 1)
            Log<MultipleReactantInEedfCol>::Warning(*this);
        if (isReverse && rhsStates.size() > 1)
            Log<MultipleProductsInReverse>::Error(*this);
    }

    EedfCollision::~EedfCollision()
    {
    }

    /** \todo What is intended here? Should we not just compare the
     *        lists of participants and the type below (after ordering),
     *        as well as the type? It may be better to store the products
     *        (and the LHS) as a map<name,stoich_coefficients> then. Then
     *        the below could be implemented as:
     *
     *        type()=o.type() && lhs()=o.lhs() && rhs==o.rhs()
     *
     */
    bool EedfCollision::is_same_as(CollisionType other_type, const EedfState* other_target, const StateVector& other_products) const
    {
        if (type != other_type) return false;
        if (getTarget() != other_target) return false;

        // Comparing pointers since there is only one instantiation of each state.
        for (const auto *state : other_products)
        {
            if (std::find(m_rhsHeavyStates.begin(), m_rhsHeavyStates.end(), state) == m_rhsHeavyStates.end())
                return false;
        }

        return true;
    }

    const EedfCollision::EedfState *EedfCollision::getTarget() const
    {
        return m_lhsHeavyStates.front();
    }
    EedfCollision::EedfState *EedfCollision::getTarget()
    {
        return m_lhsHeavyStates.front();
    }

    std::ostream &operator<<(std::ostream &os, const EedfCollision &collision)
    {
        os << "e + " << *collision.getTarget()
           << (collision.isReverse ? " <->" : " ->");

        if (collision.type != CollisionType::attachment) os << " e +";
        if (collision.type == CollisionType::ionization) os << " e +";

        for (uint32_t i = 0; i < collision.m_rhsHeavyStates.size(); ++i)
        {
            os << ' ' << *collision.m_rhsHeavyStates[i];
            if (i < collision.m_rhsHeavyStates.size() - 1) os << " +";
        }
        os << ", " << collision.typeAsString();

        return os;
    }

    void EedfCollision::superElastic(const Vector &energyData, Vector &result) const
    {
        if (!isReverse)
        {
            Log<SuperElasticForNonReverse>::Error(*this);
        }

        result.resize(energyData.size());

        Vector superElasticEnergies = energyData.array() + crossSection->threshold;

        crossSection->interpolate(superElasticEnergies, result);

        if (getTarget()->statisticalWeight <= 0.) Log<NoStatWeight>::Error(*getTarget());
        if (m_rhsHeavyStates[0]->statisticalWeight <= 0.) Log<NoStatWeight>::Error(*m_rhsHeavyStates[0]);

        const double swRatio = getTarget()->statisticalWeight / m_rhsHeavyStates[0]->statisticalWeight;

        uint32_t minIndex = 0;

        if (energyData[0] == 0)
            ++minIndex;

        for (uint32_t i = minIndex; i < result.size(); ++i)
        {
            result[i] *= swRatio * (1 + (crossSection->threshold / energyData[i]));
        }
        if (energyData[0] == 0) result[0] = 0;
    }

    CollPower EedfCollision::evaluateConservativePower(const Vector &eedf)
    {
        const Grid *grid = crossSection->getGrid();
        const uint32_t n = grid->cellNumber;

        CollPower collPower;

        Vector cellCrossSection(n);

        for (uint32_t i = 0; i < n; ++i)
        {
            cellCrossSection[i] = .5 * ((*crossSection)[i] + (*crossSection)[i + 1]);
        }

        auto lmin = static_cast<uint32_t>(crossSection->threshold / grid->step);

        const double factor = sqrt(2 * Constant::electronCharge / Constant::electronMass);

        double ineSum = 0;

        for (uint32_t i = lmin; i < n; ++i)
        {
            ineSum += eedf[i] * grid->getCell(i) * cellCrossSection[i];
        }

        collPower.ine = -factor * getTarget()->density * grid->step * grid->getNode(lmin) * ineSum;

        if (isReverse)
        {
            const double statWeightRatio = getTarget()->statisticalWeight / m_rhsHeavyStates[0]->statisticalWeight;

            double supSum = 0;

            for (uint32_t i = lmin; i < n; ++i)
            {
                supSum += eedf[i - lmin] * grid->getCell(i) * cellCrossSection[i];
            }

            collPower.sup +=
                    factor * statWeightRatio * m_rhsHeavyStates[0]->density * grid->step * grid->getNode(lmin) * supSum;
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

        for (uint32_t i = 0; i < n; ++i)
        {
            cellCrossSection[i] = .5 * ((*crossSection)[i] + (*crossSection)[i + 1]);
        }

        auto lmin = static_cast<uint32_t>(crossSection->threshold / grid->step);

        if (type == CollisionType::ionization)
        {

            if (ionizationOperatorType == IonizationOperatorType::equalSharing)
            {
                double sumOne = 0., sumTwo = 0., sumThree = 0.;

                for (uint32_t i = lmin - 1; i < n; ++i)
                {
                    sumOne += grid->getCell(i) * grid->getCell(i) * cellCrossSection[i] * eedf[i];
                }

                for (uint32_t i = 1 + lmin; i < n; i += 2)
                {
                    const double term = grid->getCell(i) * cellCrossSection[i] * eedf[i];

                    sumTwo += term;
                    sumThree += grid->getCell(i) * term;
                }

                collPower.ine = -factor * getTarget()->density * grid->step *
                                (sumOne + 2 * grid->getCell(lmin) * sumTwo - 2 * sumThree);

            }
            else if (ionizationOperatorType == IonizationOperatorType::oneTakesAll)
            {
                double sum = 0.;

                for (uint32_t i = lmin - 1; i < n; ++i) {
                    sum += grid->getCell(i) * cellCrossSection[i] * eedf[i];
                }

                collPower.ine = -factor * getTarget()->density * grid->step * grid->getCell(lmin - 1) * sum;
            }
            else if (ionizationOperatorType == IonizationOperatorType::sdcs)
            {
                double w = OPBParameter;

                if (w == 0) w = crossSection->threshold;

                Vector TICS = Vector::Zero(grid->cellNumber);
                Vector auxVec = 1. / (1. + (grid->getCells().cwiseAbs2().array() / (w * w)));

                for (int64_t k = 1; k < n; ++k) {
                    auxVec[k] = auxVec[k] + auxVec[k - 1];
                    int64_t kmax = (k + 1 - lmin) / 2;

                    if (kmax > 0)
                        TICS[k] += cellCrossSection[k] * auxVec[kmax-1] /
                                   (w * atan((grid->getCell(static_cast<uint32_t>(k)) - crossSection->threshold) / (2 * w)));
                }

                collPower.ine = -factor * getTarget()->density * grid->getCell(lmin) * grid->step *
                                eedf.cwiseProduct(grid->getCells().cwiseProduct(grid->step * TICS)).sum();
            }
        }
        else if (type == CollisionType::attachment)
        {
            double sum = 0.;

            for (uint32_t i = lmin; i < n; ++i)
            {
                sum += eedf[i] * grid->getCell(i) * grid->getCell(i) * cellCrossSection[i];
            }

            collPower.ine = -factor * getTarget()->density * grid->step * sum;
        }
        /// \todo For other Collision types, collPower is unitialized at this point

        return collPower;
    }

    RateCoefficient EedfCollision::evaluateRateCoefficient(const Vector &eedf)
    {
        const double factor = std::sqrt(2. * Constant::electronCharge / Constant::electronMass);
        const Grid *grid = crossSection->getGrid();

        const uint32_t nNodes = grid->cellNumber + 1;
        const uint32_t nCells = grid->cellNumber;

        const auto lmin = static_cast<uint32_t>(crossSection->threshold / grid->step);

        const Vector cellCrossSection =
                .5 * (crossSection->segment(lmin, nNodes - 1 - lmin) + crossSection->tail(nNodes - 1 - lmin));

        ineRateCoeff = factor * grid->step * cellCrossSection.cwiseProduct(grid->getCells().tail(nCells - lmin)).dot(
                eedf.tail(nCells - lmin));

        if (isReverse)
        {
            const double tStatWeight = getTarget()->statisticalWeight;
            const double pStatWeight = m_rhsHeavyStates[0]->statisticalWeight;

            if (tStatWeight <= 0.) Log<NoStatWeight>::Error(*getTarget());
            if (pStatWeight <= 0.) Log<NoStatWeight>::Error(*m_rhsHeavyStates[0]);

            const double statWeightRatio = tStatWeight / pStatWeight;

            supRateCoeff = factor * statWeightRatio * grid->step *
                          cellCrossSection.cwiseProduct(grid->getCells().tail(nCells - lmin)).dot(
                                  eedf.head(nCells - lmin));
        }

        return {this, ineRateCoeff, supRateCoeff};
    }

    std::string EedfCollision::typeAsString() const
    {
        switch(type)
        {
            case CollisionType::effective:
                return "Effective";
            case CollisionType::elastic:
                return "Elastic";
            case CollisionType::excitation:
                return "Excitation";
            case CollisionType::vibrational:
                return "Vibrational";
            case CollisionType::rotational:
                return "Rotational";
            case CollisionType::ionization:
                return "Ionization";
            case CollisionType::attachment:
                return "Attachment";
            default:
                return "";
        }
    }

} // namespace loki
