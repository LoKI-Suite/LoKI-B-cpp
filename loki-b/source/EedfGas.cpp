//
// Created by daan on 2-5-19.
//

#include "LoKI-B/EedfGas.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/EedfCollision.h"
#include "LoKI-B/Log.h"

namespace loki
{

EedfGas::EedfGas(const std::string &name)
    : GasBase(name), m_collisions(static_cast<uint8_t>(CollisionType::size)),
      m_collisionsExtra(static_cast<uint8_t>(CollisionType::size)),
      m_OPBParameter(0.)
{
}

void EedfGas::addCollision(EedfCollision *collision, bool isExtra)
{
    // add to the state's list
    EedfState *target = collision->getTarget();
    (isExtra ? m_state_collisionsExtra[target] : m_state_collisions[target]).emplace_back(collision);

    // add to the gas' list
    (isExtra ? this->m_collisionsExtra[static_cast<uint8_t>(collision->type())]
             : this->m_collisions[static_cast<uint8_t>(collision->type())])
        .emplace_back(collision);
}

EedfGas::~EedfGas()
{
}

CrossSection *EedfGas::elasticCrossSectionFromEffective(Grid *energyGrid)
{
    if (collisions(CollisionType::effective).empty())
        Log<Message>::Error("Could not find effective cross section for gas " + name + ".");

    EedfCollision *eff = collisions(CollisionType::effective)[0].get();
    const Vector &rawEnergies = eff->crossSection->lookupTable().x();
    const Vector &rawEff = eff->crossSection->lookupTable().y();

    Vector rawEl = rawEff; // copy raw effective into raw elastic

    if (m_effectivePopulations.empty())
    {
        m_effectivePopulations.emplace(eff->getTarget(), 1.);
        setDefaultEffPop(eff->getTarget());
    }

    for (const auto &pair : m_effectivePopulations)
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

                if (m_effectivePopulations.count(collision->m_rhsHeavyStates[0]) == 1)
                    rawEl -= crossSection * m_effectivePopulations[collision->m_rhsHeavyStates[0]];
            }
        }
    }

    bool warn = false;

    for (uint32_t i = 0; i < rawEl.size(); ++i)
    {
        if (rawEl[i] < 0.)
        {
            rawEl[i] = 0.;
            warn = true;
        }
    }

    if (warn)
        Log<NegativeElastic>::Warning(name);

    return new CrossSection(0., energyGrid, true, rawEnergies, rawEl);
}

void EedfGas::setDefaultEffPop(EedfState *ground)
{
    // ele ground to 1
    // vib children of ele ground to Boltzmann at 300K
    // rot children of vib ground to Boltzmann at 300K

    if (!ground->children().empty())
    {
        double norm = 0;
        EedfState *childGround = ground->children()[0];

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
                childGround = child;
            }
            double effPop = child->statisticalWeight * std::exp(-child->energy / (Constant::kBeV * 300));

            m_effectivePopulations.emplace(child, effPop);
            norm += effPop;
        }
        for (auto *child : ground->children())
        {
            m_effectivePopulations[child] /= (norm / m_effectivePopulations[ground]);
        }

        setDefaultEffPop(childGround);
    }
}

std::vector<EedfGas::EedfState *> EedfGas::findStatesToUpdate()
{
    std::vector<EedfState *> statesToUpdate;

    for (auto *chargeState : get_root().children())
    {
        for (auto *eleState : chargeState->children())
        {
            if (eleState->population > 0)
            {
                auto &colls = m_state_collisions[eleState];
                auto it = find_if(colls.begin(), colls.end(),
                                  [](EedfCollision *collision) { return collision->type() == CollisionType::elastic; });

                if (it == colls.end())
                    statesToUpdate.emplace_back(eleState);
            }
        }
    }

    return statesToUpdate;
}

void EedfGas::checkElasticCollisions(State *electron, Grid *energyGrid)
{
    if (isDummy())
        return;

    std::vector<EedfState *> statesToUpdate = findStatesToUpdate();

    if (!statesToUpdate.empty())
    {
        CrossSection *elasticCS = elasticCrossSectionFromEffective(energyGrid);
        std::vector<uint16_t> stoiCoeff{1, 1};

        for (auto *state : statesToUpdate)
        {
            std::vector stateVector{electron, state};
            auto *collision =
                new EedfCollision(CollisionType::elastic, stateVector, stoiCoeff, stateVector, stoiCoeff, false);

            collision->crossSection.reset(elasticCS);

            this->addCollision(collision, false);
        }
    }
}

bool EedfGas::isDummy() const
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

const GasPower& EedfGas::evaluatePower(const IonizationOperatorType ionType, const Vector &eedf)
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

const GasPower& EedfGas::getPower() const
{
    return m_power;
}

PowerTerm EedfGas::evaluateConservativePower(const CollisionVector &collisionVector, const Vector &eedf) const
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

PowerTerm EedfGas::evaluateNonConservativePower(const CollisionVector &collisionVector,
                                                const IonizationOperatorType ionType, const Vector &eedf) const
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

} // namespace loki
