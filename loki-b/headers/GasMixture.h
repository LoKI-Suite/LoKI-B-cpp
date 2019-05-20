//
// Created by daan on 14-5-19.
//

#ifndef LOKI_CPP_GASMIXTURE_H
#define LOKI_CPP_GASMIXTURE_H

// DONE: Fill the Gas, State and Collision structures.
// DONE: Implement the templated approach to gasses and states
// TODO: Think of how we can parse the LXCat data. Especially for states, do we need an
//  auxiliary structure like StateEntry or can we immediately create the state?
//  If not then alter the GasMixture class (replace StateEntry by State).
// TODO: Implement LXCat file parsing.
// TODO: Implement Property file / function parsing.

#include "EedfGas.h"
#include "EedfState.h"
#include "Collision.h"
#include "Setup.h"
#include "Traits.h"

#include <vector>

namespace loki {
    template <typename TraitType>
    struct GasMixture {
        std::vector<typename Trait<TraitType>::Gas> gasses;

        explicit GasMixture(const ElectronKineticsSetup &setup) {
            // use 'setup.LXCatFiles' to load collisions and (empty)
            // gasses and states
            // use 'setup.LXCatFilesExtra' to load extra collisions
            // use 'setup.gasProperties' to fill the gas structures
            // use 'setup.stateProperties' to fill the state structures
        }

        typename Trait<TraitType>::Gas &
        addGas(const std::string &name) {
            auto it = std::find(gasses.begin(), gasses.end(), name);

            if (it == gasses.end()) {
                return gasses.emplace_back(name);
            }

            return *it;
        }

        typename Trait<TraitType>::State &
        addState(const StateEntry &entry) {
            auto &gas = addGas(entry.gasName);

            auto *state = addState(&gas, entry);

            for (uint8_t lvl = electronic; lvl < entry.level; ++lvl) {
                state = addState(state, entry);
            }

            return *state;
        }

        void print() {
            for (int i = 0; i < gasses.size(); ++i) {
                gasses[i].print();
            }
        }

    private:
        typename Trait<TraitType>::State *
        addState(typename Trait<TraitType>::Gas * gas, const StateEntry &entry) {

            auto it = std::find_if(gas->stateTree.begin(), gas->stateTree.end(),
                                   [&entry](typename Trait<TraitType>::State *state){return state->isAncestorOf(entry);});

            if (it == gas->stateTree.end()) {
                auto *state = gas->states.emplace_back(new typename Trait<TraitType>::State(entry, gas));
                return gas->stateTree.emplace_back(state);
            }

            return *it;
        }

        typename Trait<TraitType>::State *
        addState(typename Trait<TraitType>::State *parent, const StateEntry &entry) {

            auto it = std::find_if(parent->children.begin(), parent->children.end(),
                                   [&entry](typename Trait<TraitType>::State *child){return child->isAncestorOf(entry);});

            if (it == parent->children.end()) {
                auto *state = parent->gas->states.emplace_back(new typename Trait<TraitType>::State(entry, parent->gas, parent));
                return parent->children.emplace_back(state);
            }

            return *it;
        }
    };
}


#endif //LOKI_CPP_GASMIXTURE_H
