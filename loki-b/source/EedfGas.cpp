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

    void EedfGas::checkElasticCollisions(Grid *energyGrid) {
        if (collisions.empty())
            return;

        std::vector<EedfState *> statesToUpdate;

        for (auto *eleState : stateTree) {
            if (eleState->population > 0) {
                auto it = find_if(eleState->collisions.begin(), eleState->collisions.end(),
                                  [](EedfCollision *collision) {
                                      return collision->type == CollisionType::elastic;
                                  });

                if (it == eleState->collisions.end())
                    statesToUpdate.emplace_back(eleState);
            }
        }

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

    bool EedfGas::isDummy() {
        for (const auto &vec : collisions)
            if (!vec.empty()) return false;

        return true;
    }
}