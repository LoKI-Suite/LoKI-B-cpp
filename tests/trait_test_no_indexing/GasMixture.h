//
// Created by daan on 16-5-19.
//

#ifndef LOKI_CPP_GASMIXTURE_H
#define LOKI_CPP_GASMIXTURE_H

#include "Basic.h"
#include "Traits.h"
#include "Gas.h"
#include "State.h"

#include <string>
#include <vector>
#include <algorithm>

template <typename TraitType>
struct GasMixture {
    std::vector<typename Trait<TraitType>::Gas> gasses;

    typename Trait<TraitType>::Gas &addGas(const std::string &name) {
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

    typename Trait<TraitType>::State *addState(typename Trait<TraitType>::Gas * gas, const StateEntry &entry) {

        auto it = std::find_if(gas->stateTree.begin(), gas->stateTree.end(),
                [&entry](typename Trait<TraitType>::State *state){return state->isAncestorOf(entry);});

        if (it == gas->stateTree.end()) {
            auto *state = gas->states.emplace_back(new typename Trait<TraitType>::State(entry, gas));
            return gas->stateTree.emplace_back(state);
        }

        return *it;
    }

    typename Trait<TraitType>::State *addState(typename Trait<TraitType>::State *parent, const StateEntry &entry) {

        auto it = std::find_if(parent->children.begin(), parent->children.end(),
                               [&entry](typename Trait<TraitType>::State *child){return child->isAncestorOf(entry);});

        if (it == parent->children.end()) {
            auto *state = parent->gas->states.emplace_back(new typename Trait<TraitType>::State(entry, parent->gas, parent));
            return parent->children.emplace_back(state);
        }

        return *it;
    }
};

#endif //LOKI_CPP_GASMIXTURE_H
