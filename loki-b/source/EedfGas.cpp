//
// Created by daan on 2-5-19.
//

#include "EedfGas.h"
#include "Constant.h"

namespace loki {
    EedfGas::EedfGas(const std::string &name) : Gas(name), collisions((uint8_t) Enumeration::CollisionType::size),
                                                extraCollisions((uint8_t) Enumeration::CollisionType::size) {}

    void EedfGas::addCollision(EedfCollision *collision, bool isExtra) {

        (isExtra ? this->extraCollisions[(uint8_t) collision->type] :
         this->collisions[(uint8_t) collision->type]).emplace_back(collision);
    }

    EedfGas::~EedfGas() {
        for (const auto &vec : collisions)
            for (auto *collision : vec)
                delete collision;

        for (const auto &vec : extraCollisions)
            for (auto *collision : vec)
                delete collision;
    }

    CrossSection *EedfGas::elasticCrossSectionFromEffective(Grid *energyGrid) {
        if (collisions[(uint8_t) CollisionType::effective].empty())
            Log<Message>::Error("Could not find effective cross section for gas " + name + ".");

        EedfCollision *eff = collisions[(uint8_t) CollisionType::effective][0];
        Vector &rawEff = eff->crossSection->raw();
        Vector &rawEnergies = eff->crossSection->energies();

        Vector rawEl = rawEff; // copy raw effective into raw elastic

        if (effectivePopulations.empty()) {
            effectivePopulations.emplace(eff->target, 1.);
            setDefaultEffPop(eff->target);
        }

        for (const auto &pair : effectivePopulations) {
            for (const auto *collision : pair.first->collisions) {
                if (collision->type == CollisionType::effective)
                    continue;

                Vector crossSection;
                collision->crossSection->interpolate(rawEnergies, crossSection);

                rawEl -= crossSection * pair.second;

                if (collision->isReverse) {
                    collision->superElastic(rawEnergies, crossSection);

                    if (effectivePopulations.count(collision->products[0]) == 1)
                        rawEl -= crossSection * effectivePopulations[collision->products[0]];
                }
            }
        }

        bool warn = false;

        for (uint32_t i = 0; i < rawEl.size(); ++i) {
            if (rawEl[i] < 0.) {
                rawEl[i] = 0.;
                warn = true;
            }
        }

        if (warn) Log<NegativeElastic>::Warning(name);

        return new CrossSection(0., energyGrid, rawEnergies, rawEl);
    }

    void EedfGas::setDefaultEffPop(EedfState *ground) {
        // ele ground to 1
        // vib children of ele ground to Boltzmann at 300K
        // rot children of vib ground to Boltzmann at 300K

        if (!ground->children.empty()) {
            double norm = 0;
            EedfState *childGround = ground->children[0];

            for (auto *child : ground->children) {
                if (child->energy < 0) {
                    Log<NoEnergy>::Error(*child);
                } else if (child->statisticalWeight < 0) {
                    Log<NoStatWeight>::Error(*child);
                } else if (child->energy < childGround->energy) {
                    childGround = child;
                }
                double effPop = child->statisticalWeight * std::exp(-child->energy / (Constant::kBeV * 300));

                effectivePopulations.emplace(child, effPop);
                norm += effPop;
            }
            for (auto *child : ground->children) {
                effectivePopulations[child] /= (norm / effectivePopulations[ground]);
            }

            setDefaultEffPop(childGround);
        }
    }

    void EedfGas::findStatesToUpdate(const std::vector<EedfState *> &stateStructure,
                                     std::vector<EedfState *> &statesToUpdate) {
        for (auto *eleState : stateStructure) {
            if (eleState->population > 0) {
                auto it = find_if(eleState->collisions.begin(), eleState->collisions.end(),
                                  [](EedfCollision *collision) {
                                      return collision->type == CollisionType::elastic;
                                  });

                if (it == eleState->collisions.end())
                    statesToUpdate.emplace_back(eleState);
            }
        }
    }

    void EedfGas::checkElasticCollisions(Grid *energyGrid) {
        if (isDummy()) return;

        std::vector<EedfState *> statesToUpdate;

        findStatesToUpdate(stateTree, statesToUpdate);
        findStatesToUpdate(ionicStates, statesToUpdate);

        if (!statesToUpdate.empty()) {
            CrossSection *elasticCS = elasticCrossSectionFromEffective(energyGrid);
            std::vector<uint16_t> stoiCoeff{1};

            for (auto *state : statesToUpdate) {
                std::vector stateVector{state};
                auto *collision = new EedfCollision(CollisionType::elastic, stateVector, stateVector,
                                                    stoiCoeff, false);

                collision->crossSection = elasticCS;

                state->addCollision(collision, false);
                this->addCollision(collision, false);
            }
        }
    }

    void EedfGas::checkCARConditions() {
        if (electricQuadrupoleMoment < 0)
            Log<NoElectricQuadMoment>::Error(name);

        if (rotationalConstant < 0)
            Log<NoRotationalConstant>::Error(name);

        if (!collisions[(uint8_t) CollisionType::rotational].empty())
            Log<RotCollisionInCARGas>::Error(name);
    }

    bool EedfGas::isDummy() {
        for (const auto &vec : collisions)
            if (!vec.empty()) return false;

        return true;
    }

    const GasPower &EedfGas::getPower() {
        return power;
    }

    void EedfGas::evaluatePower(const IonizationOperatorType ionType, const Vector &eedf) {

        CollPower collisionPower;

        auto *powerPtr = (double *) &power;

        for (auto i = (uint8_t) CollisionType::excitation; i <= (uint8_t) CollisionType::attachment; ++i) {

            if (i == (uint8_t) CollisionType::ionization) {
                if (ionType == IonizationOperatorType::conservative) {
                    collisionPower = evaluateConservativePower(collisions[i], eedf);
                } else {
                    collisionPower = evaluateNonConservativePower(collisions[i], ionType, eedf);
                }

                uint8_t index = 9;

                powerPtr[index] = collisionPower.ine;
            } else if (i == (uint8_t) CollisionType::attachment) {
                collisionPower = evaluateNonConservativePower(collisions[i], ionType, eedf);

                powerPtr[10] = collisionPower.ine;
            } else {
                collisionPower = evaluateConservativePower(collisions[i], eedf);

                uint8_t baseIndex = (i - (uint8_t) CollisionType::excitation) * 3;

                powerPtr[baseIndex] = collisionPower.ine;
                powerPtr[baseIndex + 1] = collisionPower.sup;
                powerPtr[baseIndex + 2] = collisionPower.ine + collisionPower.sup;
            }
        }
    }

    CollPower EedfGas::evaluateConservativePower(std::vector<EedfCollision *> &collisionVector, const Vector &eedf) {
        CollPower collPower;

        for (auto *collision : collisionVector) {
            collPower += collision->evaluateConservativePower(eedf);
        }

        return collPower;
    }

    CollPower EedfGas::evaluateNonConservativePower(std::vector<EedfCollision *> &collisionVector,
                                                    const IonizationOperatorType ionType, const Vector &eedf) {
        CollPower collPower;

        for (auto *collision : collisionVector) {
            collPower += collision->evaluateNonConservativePower(eedf, ionType, OPBParameter);
        }

        return collPower;
    }

}