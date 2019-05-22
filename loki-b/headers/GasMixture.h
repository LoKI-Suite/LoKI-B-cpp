//
// Created by daan on 14-5-19.
//

#ifndef LOKI_CPP_GASMIXTURE_H
#define LOKI_CPP_GASMIXTURE_H

// DONE: Fill the Gas, State and Collision structures.
// DONE: Implement the templated approach to gasses and states
// DONE: Think of how we can parse the LXCat data. Especially for states, do we need an
//  auxiliary structure like StateEntry or can we immediately create the state?
//  If not then alter the GasMixture class (replace StateEntry by State).
// DONE: Implement LXCat file parsing.
// TODO: Implement Property file / function parsing.

#include "EedfGas.h"
#include "EedfState.h"
#include "EedfCollision.h"
#include "Collision.h"
#include "Setup.h"
#include "Traits.h"
#include "Log.h"

#include <vector>
#include <set>
#include <fstream>
#include <regex>

// TODO: link ionic states

namespace loki {
    // TODO: move these structures
    struct RawLXCatEntry {
        std::string reactants,
                isReverse,
                products,
                type,
                threshold;
    };

    struct CollisionEntry {
        std::vector<StateEntry> reactants, products;
        std::vector<uint16_t> stoiCoeff;
        Enumeration::CollisionType type;
        double threshold;
        bool isReverse;
//        std::vector<std::pair<double, double>> rawCrossSection;
    };

    template<typename TraitType>
    struct GasMixture {
        std::vector<typename Trait<TraitType>::Gas> gasses;

        GasMixture() = default;

        virtual ~GasMixture() = default;

        void print() {
            for (int i = 0; i < gasses.size(); ++i) {
                gasses[i].print();
            }
        }

    protected:

        typename Trait<TraitType>::Gas &
        addGas(const std::string &name) {
            auto it = std::find(gasses.begin(), gasses.end(), name);

            if (it == gasses.end()) {
                return gasses.emplace_back(name);
            }

            return *it;
        }

        typename Trait<TraitType>::State *
        addState(const StateEntry &entry) {
            auto &gas = addGas(entry.gasName);

            auto *state = addState(&gas, entry);

            for (uint8_t lvl = electronic; lvl < entry.level; ++lvl) {
                state = addState(state, entry);
            }

            return state;
        }

        typename Trait<TraitType>::State *
        addState(typename Trait<TraitType>::Gas *gas, const StateEntry &entry) {

            auto it = std::find_if(gas->stateTree.begin(), gas->stateTree.end(),
                                   [&entry](typename Trait<TraitType>::State *state) {
                                       return *state >= entry;
                                   });

            if (it == gas->stateTree.end()) {
                auto *state = gas->states.emplace_back(new typename Trait<TraitType>::State(entry, gas));
                return gas->stateTree.emplace_back(state);
            }

            return *it;
        }

        typename Trait<TraitType>::State *
        addState(typename Trait<TraitType>::State *parent, const StateEntry &entry) {

            auto it = std::find_if(parent->children.begin(), parent->children.end(),
                                   [&entry](typename Trait<TraitType>::State *child) {
                                       return *child >= entry;
                                   });

            if (it == parent->children.end()) {
                auto *state = parent->gas->states.emplace_back(
                        new typename Trait<TraitType>::State(entry, parent->gas, parent));
                return parent->children.emplace_back(state);
            }

            return *it;
        }


        typename Trait<TraitType>::Collision *createCollision(CollisionEntry &entry) {
            std::vector<typename Trait<TraitType>::State *> reactants, products;
            std::set<typename Trait<TraitType>::Gas *> targetGasses;

            for (auto &stateEntry : entry.reactants)
            {
                auto *state = reactants.emplace_back(addState(stateEntry));
                targetGasses.insert(state->gas);
            }

            if (targetGasses.size() != 1)
                Log<Message>::Error("Multiple target gasses in a single collision.");

            for (auto &stateEntry : entry.products)
            {
                products.emplace_back(addState(stateEntry));
            }

            return new typename Trait<TraitType>::Collision(entry.type, reactants, products,
                                                            entry.stoiCoeff, entry.isReverse);
        }
    };
}


#endif //LOKI_CPP_GASMIXTURE_H
