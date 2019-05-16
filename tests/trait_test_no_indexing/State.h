//
// Created by daan on 16-5-19.
//

#ifndef LOKI_CPP_STATE_H
#define LOKI_CPP_STATE_H

#include "Basic.h"
#include "Traits.h"

#include <cstdint>
#include <string>
#include <vector>
#include <iostream>

template<typename TraitType>
struct State {
    typename Trait<TraitType>::State *parent{nullptr};
    typename Trait<TraitType>::Gas *gas{nullptr};

    StateLevel level;
    std::string e;
    uint16_t v, J;
    int16_t charge;

    std::vector<typename Trait<TraitType>::State *> children;

    State(const StateEntry &entry, typename Trait<TraitType>::Gas *gas, typename Trait<TraitType>::State * parent = nullptr)
            : e(entry.e), v(entry.v), J(entry.J), charge(entry.charge), parent(parent),
              gas(gas) {

        if (parent == nullptr)
            level = electronic;
        else
            level = (StateLevel)(parent->level + 1);
    }

    virtual ~State() = default;

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

    bool isAncestorOf(const typename Trait<TraitType>::State &other) {
        if (level > other.level || charge != other.charge) return false;

        if (e == other.e) {
            if (level == electronic) return true;
        } else {
            return false;
        }

        if (v == other.v) {
            if (level == vibrational) return true;
        } else {
            return false;
        }

        if (J == other.J) {
            if (level == rotational) return true;
        } else {
            return false;
        }

        return true;
    }

    void print() const {
        std::cout << "charge = " << charge << " e = " << e;
        if (level >= vibrational) std::cout << " v = " << v;
        if (level >= rotational) std::cout << " J = " << J;
        std::cout << std::endl;

        for (auto state : children) {
            if (state != nullptr) state->print();
        }
    }
};

struct EedfState : public State<Boltzmann> {
    ChemState *chemEquivalent{nullptr};

    void someExtraFunction() {
        std::cout << "I am an EedfState." << std::endl;
    }

    explicit EedfState(const StateEntry &entry, EedfGas *gas = nullptr, EedfState *parent = nullptr)
            : State<Boltzmann>(entry, gas, parent) {}
};

struct ChemState : public State<Chemistry> {
    EedfState *eedfEquivalent{nullptr};

    void someExtraFunction() {
        std::cout << "I am a ChemState." << std::endl;
    }

    explicit ChemState(const StateEntry &entry, ChemGas *gas, ChemState *parent = nullptr)
            : State<Chemistry>(entry, gas, parent) {}
};

#endif //LOKI_CPP_STATE_H
