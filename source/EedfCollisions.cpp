/** \file
 *
 *  Implementation of LoKI-B's code for electron-impact collisions and
 *  containers (per gas and for a mixture) of such collisions.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2020 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \date   21 May 2019
 */

#include "LoKI-B/EedfCollisions.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/GridOps.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/StateEntry.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>

namespace loki
{

namespace
{

template <class VectorType>
VectorType remove_electron_entries(const Collision::StateVector &parts, const VectorType &vec)
{
    assert(parts.size() == vec.size());
    VectorType result;
    // now remove all electrons from the right-hand side (since EedfCollision wants that).
    // beware of iterator invalidation
    for (unsigned i = 0; i != parts.size(); ++i)
    {
        if (parts[i]->gas().name() != "e")
        {
            result.push_back(vec[i]);
        }
    }
    return result;
}

} // namespace

EedfCollision::EedfCollision(CollisionType type, const StateVector &lhsStates, const CoeffVector &lhsCoeffs,
                             const StateVector &rhsStates, const CoeffVector &rhsCoeffs, bool isReverse)
    : Collision(type, lhsStates, lhsCoeffs, rhsStates, rhsCoeffs, isReverse),
      m_lhsHeavyStates(remove_electron_entries(lhsStates, lhsStates)),
      m_lhsHeavyCoeffs(remove_electron_entries(lhsStates, lhsCoeffs)),
      m_rhsHeavyStates(remove_electron_entries(rhsStates, rhsStates)),
      m_rhsHeavyCoeffs(remove_electron_entries(rhsStates, rhsCoeffs)), m_ineRateCoeff{0.0}, m_supRateCoeff{0.0}
{
    // for the left hand side we expect 'e + X' or 'X + e'.
    assert(lhsStates.size() == lhsCoeffs.size());
    if (lhsStates.size() != 2 || lhsCoeffs[0] != 1 || lhsCoeffs[1] != 1 ||
        (lhsStates[0]->gas().name() == "e" && lhsStates[1]->gas().name() == "e") ||
        (lhsStates[0]->gas().name() != "e" && lhsStates[1]->gas().name() != "e"))
    {
        Log<Message>::Error("Expected a binary electron impact process.");
    }
    if (m_lhsHeavyCoeffs.size() != 1 || m_lhsHeavyCoeffs[0] != 1)
    {
        Log<MultipleReactantInEedfCol>::Warning(*this);
    }
    if (isReverse)
    {
        if (m_rhsHeavyStates.size() != 1 || m_rhsHeavyCoeffs[0] != 1)
            Log<MultipleProductsInReverse>::Error(*this);
    }
}

EedfCollision::~EedfCollision()
{
}

const EedfCollision::EedfState *EedfCollision::getTarget() const
{
    return m_lhsHeavyStates.front();
}

std::ostream &operator<<(std::ostream &os, const EedfCollision &collision)
{
    os << "e + " << *collision.getTarget() << (collision.isReverse() ? " <->" : " ->");

    if (collision.type() != CollisionType::attachment)
        os << " e +";
    if (collision.type() == CollisionType::ionization)
        os << " e +";

    for (EedfCollision::StateVector::size_type i = 0; i != collision.m_rhsHeavyStates.size(); ++i)
    {
        os << ' ' << *collision.m_rhsHeavyStates[i];
        if (i < collision.m_rhsHeavyStates.size() - 1)
            os << " +";
    }
    os << ", " << collision.typeAsString();

    return os;
}

void EedfCollision::superElastic(const Vector &energyData, Vector &result) const
{
    if (!isReverse())
    {
        Log<SuperElasticForNonReverse>::Error(*this);
    }

    result.resize(energyData.size());

    Vector superElasticEnergies = energyData.array() + crossSection->threshold();

    crossSection->interpolate(superElasticEnergies, result);

    if (getTarget()->statisticalWeight <= 0.)
        Log<NoStatWeight>::Error(*getTarget());
    if (m_rhsHeavyStates[0]->statisticalWeight <= 0.)
        Log<NoStatWeight>::Error(*m_rhsHeavyStates[0]);

    const double swRatio = getTarget()->statisticalWeight / m_rhsHeavyStates[0]->statisticalWeight;

    uint32_t minIndex = 0;

    if (energyData[0] == 0)
        ++minIndex;

    for (uint32_t i = minIndex; i < result.size(); ++i)
    {
        result[i] *= swRatio * (1 + (crossSection->threshold() / energyData[i]));
    }
    if (energyData[0] == 0)
        result[0] = 0;
}

PowerTerm EedfCollision::evaluateConservativePower(const Vector &eedf) const
{
    const Grid *grid = crossSection->getGrid();
    const Grid::Index n = grid->nCells();

    PowerTerm collPower;

    Vector cellCrossSection(n);

    for (uint32_t i = 0; i < n; ++i)
    {
        cellCrossSection[i] = .5 * ((*crossSection)[i] + (*crossSection)[i + 1]);
    }
    int lmin;
    if (grid->isUniform())
    {
        lmin = static_cast<uint32_t>(crossSection->threshold() / grid->du());
    } else
    {
        lmin = static_cast<uint32_t>((std::upper_bound(grid->getNodes().begin(),grid->getNodes().end(), crossSection->threshold()) - grid->getNodes().begin())) - 1;
    }
    
    double ineSum = 0;

    if (grid->isUniform())
    {
        for (uint32_t i = lmin; i < n; ++i)
        {
            ineSum += eedf[i] * grid->getCell(i) * cellCrossSection[i];
        }
        collPower.forward = -SI::gamma * getTarget()->delta() * grid->du() * grid->getNode(lmin) * ineSum;
    } else
    {
        for (uint32_t i = lmin; i < n; ++i)
        {
            ineSum += eedf[i] * grid->getCell(i) * cellCrossSection[i] * grid->duCell(i);
        }
        collPower.forward = -SI::gamma * getTarget()->delta() *  grid->getNode(lmin) * ineSum;
    }
    

    if (isReverse())
    {
        const double statWeightRatio = getTarget()->statisticalWeight / m_rhsHeavyStates[0]->statisticalWeight;

        double supSum = 0;

        if (grid->isUniform())
        {
            for (uint32_t i = lmin; i < n; ++i)
            {
                supSum += eedf[i - lmin] * grid->getCell(i) * cellCrossSection[i];
            }
            collPower.backward +=
                SI::gamma * statWeightRatio * m_rhsHeavyStates[0]->delta() * grid->du() * grid->getNode(lmin) * supSum;
        } else
        {
            for (uint32_t i = lmin; i < n; ++i)
            {
                supSum += eedf[i - lmin] * grid->getCell(i) * cellCrossSection[i] * grid->duCell(i - lmin);
            }
            collPower.backward +=
                SI::gamma * statWeightRatio * m_rhsHeavyStates[0]->delta() * grid->getNode(lmin) * supSum;
        }
        
    }

    return collPower;
}

PowerTerm EedfCollision::evaluateNonConservativePower(const Vector &eedf,
                                                      const IonizationOperatorType ionizationOperatorType,
                                                      const double OPBParameter) const
{
    const Grid *grid = crossSection->getGrid();
    const Grid::Index n = grid->nCells();

    PowerTerm collPower;

    Vector cellCrossSection(n);

    for (uint32_t i = 0; i < n; ++i)
    {
        cellCrossSection[i] = .5 * ((*crossSection)[i] + (*crossSection)[i + 1]);
    }

    auto lmin = static_cast<uint32_t>(crossSection->threshold() / grid->du());

    if (type() == CollisionType::ionization)
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

            collPower.forward = -SI::gamma * getTarget()->delta() * grid->du() *
                                (sumOne + 2 * grid->getCell(lmin) * sumTwo - 2 * sumThree);
        }
        else if (ionizationOperatorType == IonizationOperatorType::oneTakesAll)
        {
            double sum = 0.;

            for (uint32_t i = lmin - 1; i < n; ++i)
            {
                sum += grid->getCell(i) * cellCrossSection[i] * eedf[i];
            }

            collPower.forward = -SI::gamma * getTarget()->delta() * grid->du() * grid->getCell(lmin - 1) * sum;
        }
        else if (ionizationOperatorType == IonizationOperatorType::sdcs)
        {
            double w = OPBParameter;

            if (w < 0)
                w = crossSection->threshold();

            Vector TICS = Vector::Zero(grid->nCells());
            Vector auxVec = 1. / (1. + (grid->getCells().cwiseAbs2().array() / (w * w)));

            for (int64_t k = 1; k < n; ++k)
            {
                auxVec[k] = auxVec[k] + auxVec[k - 1];
                int64_t kmax = (k + 1 - lmin) / 2;

                if (kmax > 0)
                    TICS[k] +=
                        cellCrossSection[k] * auxVec[kmax - 1] /
                        (w * atan((grid->getCell(static_cast<uint32_t>(k)) - crossSection->threshold()) / (2 * w)));
            }

            collPower.forward = -SI::gamma * getTarget()->delta() * grid->getCell(lmin) * grid->du() *
                                eedf.cwiseProduct(grid->getCells().cwiseProduct(grid->du() * TICS)).sum();
        }
    }
    else if (type() == CollisionType::attachment)
    {
        double sum = 0.;

        for (uint32_t i = lmin; i < n; ++i)
        {
            sum += eedf[i] * grid->getCell(i) * grid->getCell(i) * cellCrossSection[i];
        }

        collPower.forward = -SI::gamma * getTarget()->delta() * grid->du() * sum;
    }
    /// \todo For other Collision types, collPower is unitialized at this point

    return collPower;
}

RateCoefficient EedfCollision::evaluateRateCoefficient(const Vector &eedf)
{
    const Grid *grid = crossSection->getGrid();
    if (crossSection->threshold() > grid->uMax())
    {
        return {this, 0.0, 0.0};
    }

    const Grid::Index nNodes = grid->nCells() + 1;
    const Grid::Index nCells = grid->nCells();

    int lmin;
    if (grid->isUniform())
    {
        lmin = static_cast<uint32_t>(crossSection->threshold() / grid->du());
    } else
    {
        lmin = static_cast<Grid::Index>(std::upper_bound(grid->getNodes().begin(),grid->getNodes().end(), crossSection->threshold()) - grid->getNodes().begin()) - 1;
    }
    

    const Vector cellCrossSection =
        .5 * (crossSection->segment(lmin, nNodes - 1 - lmin) + crossSection->tail(nNodes - 1 - lmin));

    if (grid->isUniform())
    {
        m_ineRateCoeff = SI::gamma * grid->du() *
                     cellCrossSection.cwiseProduct(grid->getCells().tail(nCells - lmin)).dot(eedf.tail(nCells - lmin));

    } else
    {
        m_ineRateCoeff = SI::gamma * (grid->duCells().tail(nCells - lmin)).dot(cellCrossSection) * 
                     grid->getCells().tail(nCells - lmin).dot(eedf.tail(nCells - lmin));
    }

    if (isReverse())
    {
        const double tStatWeight = getTarget()->statisticalWeight;
        const double pStatWeight = m_rhsHeavyStates[0]->statisticalWeight;

        if (tStatWeight <= 0.)
            Log<NoStatWeight>::Error(*getTarget());
        if (pStatWeight <= 0.)
            Log<NoStatWeight>::Error(*m_rhsHeavyStates[0]);

        const double statWeightRatio = tStatWeight / pStatWeight;
        
        if (grid->isUniform())
        {
            m_supRateCoeff =
                SI::gamma * statWeightRatio * grid->du() *
                cellCrossSection.cwiseProduct(grid->getCells().tail(nCells - lmin)).dot(eedf.head(nCells - lmin)); 
        } else
        {
            m_supRateCoeff =
                SI::gamma * statWeightRatio *
                (grid->duCells().head(nCells - lmin)).dot(cellCrossSection) * 
                     grid->getCells().tail(nCells - lmin).dot(eedf.head(nCells - lmin)); 
        }
        
    }

    return {this, m_ineRateCoeff, m_supRateCoeff};
}

std::string EedfCollision::typeAsString() const
{
    switch (type())
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

EedfCollisionDataGas::EedfCollisionDataGas(const GasProperties& gasProps, const Gas &gas)
    : m_gas(gas), m_collisions(static_cast<uint8_t>(CollisionType::size)),
      m_collisionsExtra(static_cast<uint8_t>(CollisionType::size)),
      m_OPBParameter(gasProps.get("OPBParameter",gas.name(),-1.0,true))
{
    if (m_OPBParameter>0)
    {
        Log<Message>::Notify("Set OPB parameter for gas "
                             + m_gas.name() + " to " + std::to_string(m_OPBParameter));
    }
}

void EedfCollisionDataGas::addCollision(EedfCollision *collision, bool isExtra)
{
    const State *target = collision->getTarget();
    if (isExtra)
    {
        m_state_collisionsExtra[target].emplace_back(collision);
        m_collisionsExtra[static_cast<uint8_t>(collision->type())].emplace_back(collision);
    }
    else
    {
        m_state_collisions[target].emplace_back(collision);
        m_collisions[static_cast<uint8_t>(collision->type())].emplace_back(collision);
    }
}

void EedfCollisionDataGas::checkElasticCollisions(const State *electron, const Grid *energyGrid, const EffectivePopulationsMap& effectivePopulation)
{
    if (isDummy())
        return;

    std::vector<const State *> statesToUpdate = findStatesToUpdate();

    if (!statesToUpdate.empty())
    {
        CrossSection *elasticCS = elasticCrossSectionFromEffective(energyGrid,effectivePopulation);
        std::vector<uint16_t> stoiCoeff{1, 1};

        for (const auto *state : statesToUpdate)
        {
            Log<Message>::Notify("Effective cross section: installing elastic cross section for state " + (std::stringstream{}<<*state).str());
            std::vector stateVector{electron, state};
            auto *collision =
                new EedfCollision(CollisionType::elastic, stateVector, stoiCoeff, stateVector, stoiCoeff, false);

            collision->crossSection.reset(elasticCS);

            this->addCollision(collision, false);
        }
    }
}

CrossSection *EedfCollisionDataGas::elasticCrossSectionFromEffective(const Grid *energyGrid, const EffectivePopulationsMap& effectivePopulationsCustom)
{
    /** \todo What happens / should happen if more than one effective collision
     *  is specified? Is that an input error? The code below only uses the first.
     */
    if (collisions(CollisionType::effective).empty())
        Log<Message>::Error("Could not find effective cross section for gas " + m_gas.name() + ".");

    EedfCollision *eff = collisions(CollisionType::effective)[0].get();
    const Vector &rawEnergies = eff->crossSection->lookupTable().x();
    const Vector &rawEff = eff->crossSection->lookupTable().y();

    Vector rawEl = rawEff; // copy raw effective into raw elastic

    EffectivePopulationsMap effectivePopulationsDefault;
    if (effectivePopulationsCustom.empty())
    {
        effectivePopulationsDefault.emplace(eff->getTarget(), 1.);
        setDefaultEffPop(eff->getTarget(),effectivePopulationsDefault);
    }

    const EffectivePopulationsMap& effectivePopulations = effectivePopulationsCustom.empty() ? effectivePopulationsDefault : effectivePopulationsCustom;
    for (const auto &pair : effectivePopulations)
    {
        for (const auto &collision : m_state_collisions[pair.first])
        {
            if (collision->type() == CollisionType::effective)
                continue;

            Vector crossSection;
            collision->crossSection->interpolate(rawEnergies, crossSection);

            rawEl -= crossSection * pair.second;

            if (collision->isReverse())
            {
                collision->superElastic(rawEnergies, crossSection);

                const auto pop = effectivePopulations.find(collision->m_rhsHeavyStates[0]);
                if (pop != effectivePopulations.end())
                {
                    rawEl -= crossSection * pop->second;
                }
            }
        }
    }

    bool warn = false;

    for (Grid::Index i = 0; i < rawEl.size(); ++i)
    {
        if (rawEl[i] < 0.)
        {
            rawEl[i] = 0.;
            warn = true;
        }
    }

    if (warn)
        Log<NegativeElastic>::Warning(m_gas.name());

    return new CrossSection(0., energyGrid, true, rawEnergies, rawEl);
}

/** \todo This should probably be in Gas; this is a generic thermodynamic
 *  operation that does not depend on collisional data.
 *  \todo The 300 (K) should not be hardcoded.
 */
void EedfCollisionDataGas::setDefaultEffPop(const State *ground, EffectivePopulationsMap& effectivePopulations) const
{
    // ele ground to 1
    // vib children of ele ground to Boltzmann at 300K
    // rot children of vib ground to Boltzmann at 300K

    if (!ground->children().empty())
    {
        double norm = 0;
        State *childGround = ground->children()[0];

        /** \todo Check the algorithm. It should be possible to
         *  implement this without identifying the child ground.
         */
        for (auto *child : ground->children())
        {
            if (child->energy < 0)
            {
                Log<NoEnergy>::Error(*child);
            }
            else if (child->statisticalWeight < 0)
            {
                Log<NoStatWeight>::Error(*child);
            }
            else if (child->energy < childGround->energy)
            {
                /// \todo Explain the logic of this line. Also like this in MATLAB?
                childGround = child;
            }
            double effPop = child->statisticalWeight * std::exp(-child->energy / (Constant::kBeV * 300));
            /// \todo Explain that the ground state weight does not matter in view of the normalization below.

            effectivePopulations.emplace(child, effPop);
            norm += effPop;
        }
        for (auto *child : ground->children())
        {
            effectivePopulations[child] /= (norm / effectivePopulations[ground]);
        }

        setDefaultEffPop(childGround,effectivePopulations);
    }
}

bool EedfCollisionDataGas::isDummy() const
{
    for (const auto &vec : collisions())
    {
        if (!vec.empty())
        {
            return false;
        }
    }
    return true;
}

std::vector<const EedfCollisionDataGas::State *> EedfCollisionDataGas::findStatesToUpdate() const
{
#define NEW_FINDSTATESTOUPDATE_IMPLEMENTATION 0
#if NEW_FINDSTATESTOUPDATE_IMPLEMENTATION

    std::vector<const State *> statesToUpdate;

    /** \todo This seems wrong: an effective cross section for the neutrals will
     *  will also be used as a basis for the charged states' elastic cross section,
     *  it seems.
     */
    const auto hasElastic = [this](const loki::Gas::State *state) -> bool {
            const auto colls_it = m_state_collisions.find(state);
            if (colls_it==m_state_collisions.end())
            {
                std::stringstream ss; ss << *state;
                throw std::runtime_error("No state collisions found for state '" + ss.str() + "'.");
            }
            const auto &colls = colls_it->second;
            auto it = find_if(colls.begin(), colls.end(),
                              [](const EedfCollision *collision) { return collision->type() == CollisionType::elastic; });
            return it != colls.end();
    };

    const std::function<bool(const loki::Gas::State *)> hasElasticRecursive = [hasElastic, &hasElasticRecursive](const loki::Gas::State *state) -> bool {
        // Either the state has an elastic cross section, or all of its (grand)children with
        // nonzero population have an elastic cross section.
        if (!hasElastic(state)) {
            if (state->children().empty()) {
                return false;
            }

            // This assumes that at least one child state needs to have a nonzero population.
            for (const auto * child : state->children()) {
                if (child->population() > 0) {
                    if (!hasElasticRecursive(child)) {
                        return false;
                    }
                }
            }

            // All children have an elastic cross section.
            return true;
        } else {
            return true;
        }
    };

    for (const auto *chargeState : m_gas.get_root().children())
    {
        for (const auto *eleState : chargeState->children())
        {
            if (eleState->population() > 0 && !hasElasticRecursive(eleState))
            {
                statesToUpdate.emplace_back(eleState);
            }
        }
    }

    return statesToUpdate;

#else

    // the original implementation.

    std::vector<const State *> statesToUpdate;

    /** \todo This seems wrong: an effective cross section for the neutrals will
     *  will also be used as a basis for the charged states' elastic cross section,
     *  it seems.
     */
    for (const auto *chargeState : m_gas.get_root().children())
    {
        for (const auto *eleState : chargeState->children())
        {
            if (eleState->population() > 0)
            {
                const auto colls_it = m_state_collisions.find(eleState);
                if (colls_it==m_state_collisions.end())
                {
                    std::stringstream ss; ss << *eleState;
                    throw std::runtime_error("No state collisions found for state '" + ss.str() + "'.");
                }
                const auto &colls = colls_it->second;
                const auto it = find_if(colls.begin(), colls.end(),
                                  [](const EedfCollision *collision) { return collision->type() == CollisionType::elastic; });

                if (it == colls.end())
                    statesToUpdate.emplace_back(eleState);
            }
        }
    }

    return statesToUpdate;

#endif // NEW_FINDSTATESTOUPDATE_IMPLEMENTATION
}

const GasPower &EedfCollisionDataGas::evaluatePower(const IonizationOperatorType ionType, const Vector &eedf)
{
    m_power.ionization = (ionType == IonizationOperatorType::conservative)
                             ? evaluateConservativePower(collisions(CollisionType::ionization), eedf)
                             : evaluateNonConservativePower(collisions(CollisionType::ionization), ionType, eedf);
    m_power.attachment = evaluateNonConservativePower(collisions(CollisionType::attachment), ionType, eedf);
    m_power.excitation = evaluateConservativePower(collisions(CollisionType::excitation), eedf);
    m_power.vibrational = evaluateConservativePower(collisions(CollisionType::vibrational), eedf);
    m_power.rotational = evaluateConservativePower(collisions(CollisionType::rotational), eedf);
    return m_power;
}

const GasPower &EedfCollisionDataGas::getPower() const
{
    return m_power;
}

PowerTerm EedfCollisionDataGas::evaluateConservativePower(const CollisionVector &collisionVector,
                                                          const Vector &eedf) const
{
    PowerTerm collPower;

    for (auto &collision : collisionVector)
    {
        const Grid *grid = collision->crossSection->getGrid();

        if (collision->crossSection->threshold() > grid->uMax())
            continue;
        collPower += collision->evaluateConservativePower(eedf);
    }
    return collPower;
}

PowerTerm EedfCollisionDataGas::evaluateNonConservativePower(const CollisionVector &collisionVector,
                                                             const IonizationOperatorType ionType,
                                                             const Vector &eedf) const
{
    PowerTerm collPower;

    for (auto &collision : collisionVector)
    {
        const Grid *grid = collision->crossSection->getGrid();

        if (collision->crossSection->threshold() > grid->uMax())
            continue;

        collPower += collision->evaluateNonConservativePower(eedf, ionType, OPBParameter());
    }
    return collPower;
}

EedfCollisionDataMixture::EedfCollisionDataMixture()
{
}

EedfCollision &EedfCollisionDataMixture::addCollision(CollisionType type, const Collision::StateVector &lhsStates,
                                                      const Collision::CoeffVector &lhsCoeffs,
                                                      const Collision::StateVector &rhsStates,
                                                      const Collision::CoeffVector &rhsCoeffs, bool reverseAlso,
                                                      bool isExtra)
{
    // 1. Create the collision object (but do not configure a CrossSection
    //    object yet).
    std::unique_ptr<EedfCollision> coll_uptr{
        new Collision(type, lhsStates, lhsCoeffs, rhsStates, rhsCoeffs, reverseAlso)};
    Log<Message>::Notify(*coll_uptr);
    // 2. See if we already have a collision of the same type with the same
    //    lhs and rhs. If we do, a runtime_error is thrown.
    for (const auto &c : m_collisions)
    {
        if (c.m_coll->is_same_as(*coll_uptr))
        {
            Log<DoubleCollision>::Error(*c.m_coll);
        }
    }
    // 3. Release the collision object, store the released pointer in coll.
    //    Use coll from here on. The reason is that the ownership of the
    //    pointer is taken by the gas-specific collision container, to which
    //    it will be added below.
    Collision *coll = coll_uptr.release();
    // 4. Register the collision with the (single) gas that is the target of the
    //    binary electron-heavy process, add it to our own list of collisions,
    //    mark the process type is being present and return the (released) pointer
    //    to our caller.
    //    That the left-hand side of the process is indeed of the form 'e + X' is
    //    checked by the EedfCollision constructor, and the target 'X' is returned
    //    by its member getTarget().
    get_data_for(coll->getTarget()->gas()).addCollision(coll, isExtra);
    m_collisions.emplace_back(CollisionEntry{coll, isExtra});
    m_hasCollisions[static_cast<uint8_t>(coll->type())] = true;
    return *coll;
}

EedfCollisionDataMixture::State *EedfCollisionDataMixture::ensureState(const GasProperties& gasProps, GasMixture &composition, const StateEntry &entry)
{
    // find the state in the composition object (and let that create it if it does not
    // yet exist)
    State *state = composition.ensureState(gasProps,entry);
    // check that an entry exists in the for the state's gas in the m_data_per_gas
    // array; create such entry if that is not the case
    auto it = m_data_per_gas.begin();
    for (; it != m_data_per_gas.end(); ++it)
    {
        if (&state->gas() == &it->gas())
        {
            break;
        }
    }
    if (it == m_data_per_gas.end())
    {
        m_data_per_gas.emplace_back(gasProps,state->gas());
    }
    return state;
}

void EedfCollisionDataMixture::loadCollisionsClassic(const std::filesystem::path &file, const GasProperties& gasProps, GasMixture& composition, const Grid *energyGrid,
                                                     bool isExtra)
{
    const std::regex reParam(R"(PARAM\.:)");
    const std::regex reThreshold(R"(E = +(\S*) +eV)");
    const std::regex reProcess(R"(\[(.+?)(<->|->)(.+?), (\w+)\])");
    std::ifstream in(file);
    if (!in.is_open())
    {
        Log<FileError>::Error(file.generic_string());
        return;
    }

    std::string line;
    std::smatch mThreshold, mProcess;

    /*  Expectations: we cycle throught the file and look for a line that has "PARAM.:".
     *  When we find that, we start parsing a process description.
     *  1. Read the arguments of the parameter line, see if a threshold is specified ("E = <double> eV"),
     *     Set the threshold to 0 otherwise.
     *  2. The next line should describe the process: it must be of the form [<lhs> -> <rhs>, <type>].
     *  3. A Collision is created from the process description information
     *  4. (Maybe) a cross section object is attached to the collision object, based on subsequent lines.
     *
     */
    while (std::getline(in, line))
    {
        try
        {
            if (!std::regex_search(line, reParam))
            {
                continue;
            }
            double threshold = 0.0;
            if (std::regex_search(line, mThreshold, reThreshold))
            {
		            threshold = Parse::getValue(mThreshold[1]);
            }
            if (!std::getline(in, line) || !std::regex_search(line, mProcess, reProcess))
            {
                Log<LXCatError>::Error(file.generic_string());
            }

            try
            {
                // mProcess[1...4] represent: lhs, separator, rhs and type, and could
                // be something like:
                //  1: "He + e"
                //  2: "->" ("<->" for a two-way process)
                //  3: "He + e"
                //  4: "Elastic"

                std::vector<StateEntry> entry_lhsStates, entry_rhsStates;
                std::vector<uint16_t> entry_lhsCoeffs, entry_rhsCoeffs;

                entriesFromString(mProcess[1].str(), entry_lhsStates, &entry_lhsCoeffs);
                const bool reverseAlso = (mProcess[2].str()[0] == '<');
                entriesFromString(mProcess[3].str(), entry_rhsStates, &entry_rhsCoeffs);
                const CollisionType entry_type = getCollisionType(mProcess[4].str());

                // 1. Create vectors of pointers to the states that appear
                //    on the left and right-hand sides of the process.
                //    Create the states when necessary.
                std::vector<const State *> lhsStates;
                std::vector<const State *> rhsStates;
                for (auto &stateEntry : entry_lhsStates)
                {
                    lhsStates.emplace_back(ensureState(gasProps, composition, stateEntry));
                }
                for (auto &stateEntry : entry_rhsStates)
                {
                    rhsStates.emplace_back(ensureState(gasProps, composition, stateEntry));
                }
                Collision &coll = addCollision(entry_type, lhsStates, entry_lhsCoeffs, rhsStates, entry_rhsCoeffs,
                                               reverseAlso, isExtra);
                const bool isElasticOrEffective =
                    (coll.type() == CollisionType::effective || coll.type() == CollisionType::elastic);
                coll.crossSection.reset(new CrossSection(threshold, energyGrid, isElasticOrEffective, in));
            }
            catch (std::exception &exc)
            {
                throw std::runtime_error("Error while parsing reaction from section '" + line + "':\n" +
                                         std::string(exc.what()));
            }
        }
        catch (std::exception &exc)
        {
            throw std::runtime_error("While creating collision '" + line + "':\n" + std::string(exc.what()));
        }
    }
}

void EedfCollisionDataMixture::loadCollisionsJSON(const json_type &mcnf, const GasProperties& gasProps, GasMixture &composition, const Grid *energyGrid,
                                                  bool isExtra)
{
    Log<Message>::Notify("Started loading collisions.");

    // 1. read the states
    for (const auto &[key, value] : mcnf.at("states").items())
    {
        ensureState(gasProps, composition, entryFromJSON(key, value));
    }
    // 2. read the processes
    for (const auto &pcnf : mcnf.at("processes"))
    {
        try
        {
            const json_type &rcnf = pcnf.at("reaction");

            std::vector<uint16_t> lhsCoeffs, rhsCoeffs;

            // 1. Create vectors of pointers to the states that appear
            //    on the left and right-hand sides of the process.
            std::vector<const State *> lhsStates;
            std::vector<const State *> rhsStates;
            for (const auto &t : rcnf.at("lhs"))
            {
                const std::string &stateName(t.at("state"));
                State *state = composition.findStateById(stateName);
                if (!state)
                {
                    throw std::runtime_error("Could not find state '" + stateName + "'.");
                }
                lhsStates.push_back(state);
                lhsCoeffs.push_back(t.contains("count") ? t.at("count").get<int>() : 1);
            }
            for (const auto &t : rcnf.at("rhs"))
            {
                const std::string &stateName(t.at("state"));
                State *state = composition.findStateById(stateName);
                if (!state)
                {
                    throw std::runtime_error("Could not find state '" + stateName + "'.");
                }
                rhsStates.push_back(state);
                rhsCoeffs.push_back(t.contains("count") ? t.at("count").get<int>() : 1);
            }

            const CollisionType type = getCollisionTypeFromTypeTagArray(rcnf.at("type_tags"));
            const bool reverseAlso = rcnf.at("reversible");

            Collision &coll = addCollision(type, lhsStates, lhsCoeffs, rhsStates, rhsCoeffs, reverseAlso, isExtra);
            const bool isElasticOrEffective =
                (coll.type() == CollisionType::effective || coll.type() == CollisionType::elastic);
            coll.crossSection.reset(new CrossSection(energyGrid, isElasticOrEffective, pcnf));
        }
        catch (std::exception &exc)
        {
            throw std::runtime_error("Error while parsing reaction from section '" + pcnf.dump(1) + "':\n" +
                                     std::string(exc.what()));
        }
    }
    Log<Message>::Notify("Finished loading collisions.");
}

const EedfCollisionDataGas &EedfCollisionDataMixture::get_data_for(const Gas &gas) const
{
    for (const auto &cd : m_data_per_gas)
    {
        if (&gas == &cd.gas())
        {
            return cd;
        }
    }
    throw std::runtime_error("Could not find per-gas collision data for '" + gas.name() + "'.");
}

EedfCollisionDataGas &EedfCollisionDataMixture::get_data_for(const Gas &gas)
{
    return const_cast<EedfCollisionDataGas &>(static_cast<const EedfCollisionDataMixture &>(*this).get_data_for(gas));
}

void EedfCollisionDataMixture::evaluateTotalAndElasticCS(const Grid &grid)
{
    // 1. resize and zero-initialize the elastic and total cross sections
    m_elasticCrossSection.setZero(grid.nCells() + 1);
    m_totalCrossSection.setZero(grid.nCells() + 1);

    for (const auto &cd : data_per_gas())
    {
        if (cd.isDummy())
        {
            continue;
        }
        // 2. add elastic terms
        const double massRatio = Constant::electronMass / cd.gas().mass;
        for (auto &collision : cd.collisions(CollisionType::elastic))
        {
            m_elasticCrossSection += *collision->crossSection * (collision->getTarget()->delta() * massRatio);
            m_totalCrossSection += *collision->crossSection * collision->getTarget()->delta();
        }
        // 3. add inelastic terms (also from reverse processes, if enabled)
        //    (note: here inelastic also includes non-conserving process:
        //    ionization, attachment)
        for (auto ctype : {CollisionType::excitation, CollisionType::vibrational, CollisionType::rotational,
                           CollisionType::ionization, CollisionType::attachment})

        {
            for (auto &collision : cd.collisions(ctype))
            {
                m_totalCrossSection += *collision->crossSection * collision->getTarget()->delta();

                if (collision->isReverse())
                {
                    Vector superElastic;
                    collision->superElastic(grid.getNodes(), superElastic);
                    m_totalCrossSection += superElastic * collision->m_rhsHeavyStates[0]->delta();
                }
            }
        }
    }
    interpolateNodalToCell(grid,m_totalCrossSection,m_totalCellCrossSection);
}

void EedfCollisionDataMixture::evaluateRateCoefficients(const Grid &grid, const Vector &eedf)
{
    m_rateCoefficients.clear();
    m_rateCoefficientsExtra.clear();
    for (const auto &cd : data_per_gas())
    {
        for (auto &collVec : cd.collisions())
        {
            for (auto &collision : collVec)
            {
                if (collision->crossSection->threshold() > grid.uMax())
                {
                    continue;
                }
                m_rateCoefficients.emplace_back(collision->evaluateRateCoefficient(eedf));
            }
        }
        for (auto &collVec : cd.collisionsExtra())
        {
            for (auto &collision : collVec)
            {
                if (collision->crossSection->threshold() > grid.uMax())
                {
                    continue;
                }
                m_rateCoefficientsExtra.emplace_back(collision->evaluateRateCoefficient(eedf));
            }
        }
    }
}

} // namespace loki
