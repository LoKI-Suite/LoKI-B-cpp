//
// Created by daan on 16-5-19.
//

#ifndef LOKI_CPP_WORKING_V1_H
#define LOKI_CPP_WORKING_V1_H

#include "Basic.h"

#include <vector>
#include <algorithm>
#include <iostream>

struct Gas;

struct State {
    State *parent{nullptr};
    Gas *gas{nullptr};

    StateLevel level;
    std::string e;
    uint16_t v, J;
    int16_t charge;

    std::vector<State *> children;

    State(const StateEntry &entry, Gas *gas, State * parent = nullptr)
            : e(entry.e), v(entry.v), J(entry.J), charge(entry.charge), parent(parent),
              gas(gas) {

        if (parent == nullptr)
            level = electronic;
        else
            level = (StateLevel)(parent->level + 1);
    }

    bool isAncestorOf(const StateEntry &entry) {
        if (level > entry.level || charge != entry.charge) return false;

        if (e == entry.e) {
            if (level == electronic) return true;
        } else {
            return false;
        }

        if (v == entry.v) {
            if (level == vibrational) return true;
        } else {
            return false;
        }

        if (J == entry.J) {
            if (level == rotational) return true;
        } else {
            return false;
        }

        return true;
    }

    void print() const {
        std::cout << "level: e = " << e;
        if (level >= vibrational) std::cout << " v = " << v;
        if (level >= rotational) std::cout << " J = " << J;
        if (level == ionic) std::cout << " charge = " << charge;
        std::cout << std::endl;

        for (auto state : children) {
            state->print();
        }
    }
};

struct Gas {
    const std::string name;
    std::vector<State *> states;
    std::vector<State *> stateTree;

    bool operator==(const std::string &otherName) {
        return name == otherName;
    }

    void print() const {
        std::cout << "Gas: " << name << std::endl;
        for (auto i : stateTree) {
            i->print();
        }
    }

    explicit Gas (std::string name) : name(std::move(name)) {}

    ~Gas () {
        for (auto *state : states)
            delete state;
    }
};

struct GasMixture {
    std::vector<Gas> gasses;

    Gas &addGas(const std::string &name) {
        auto it = std::find(gasses.begin(), gasses.end(), name);

        if (it == gasses.end()) {
            return gasses.emplace_back(name);
        }

        return *it;
    }

    State &addState(const StateEntry &entry) {
        auto &gas = addGas(entry.gasName);

        auto *state = addState(&gas, entry);

        for (uint8_t lvl = electronic; lvl < entry.level; ++lvl) {
            state = addState(state, entry);
        }

        return *state;
    }

    void print() {
        for (auto & gas : gasses) {
            gas.print();
        }
    }

private:

    State *addState(State *parent, const StateEntry &entry) {
        auto &children = parent->children;
        auto &gas = parent->gas;

        auto it = std::find_if(children.begin(), children.end(),
                               [&entry](State * child){return child->isAncestorOf(entry);});

        if (it == children.end()) {
            auto *state = gas->states.emplace_back(new State(entry, gas, parent));
            return children.emplace_back(state);
        }

        return *it;
    }

    State *addState(Gas * gas, const StateEntry &entry) {
        auto it = std::find_if(gas->stateTree.begin(), gas->stateTree.end(),
                               [&entry](State * state){return state->isAncestorOf(entry);});

        if (it == gas->stateTree.end()) {
            auto *state = gas->states.emplace_back(new State(entry, gas));
            return gas->stateTree.emplace_back(state);
        }

        return *it;
    }
};

#endif //LOKI_CPP_WORKING_V1_H
