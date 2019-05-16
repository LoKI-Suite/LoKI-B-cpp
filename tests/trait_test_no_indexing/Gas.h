//
// Created by daan on 16-5-19.
//

#ifndef LOKI_CPP_GAS_H
#define LOKI_CPP_GAS_H

#include "Traits.h"
#include "Basic.h"
#include "State.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

template <typename TraitType>
struct Gas {
    const std::string name;
    std::vector<typename Trait<TraitType>::State *> states;
    std::vector<typename Trait<TraitType>::State *> stateTree;

    explicit Gas(std::string name) : name(std::move(name)) {}

    ~Gas() {
        for (auto *state : states)
            delete state;
    }

    bool operator==(const std::string &otherName) {
        return name == otherName;
    }

    void print() const {
        std::cout << "Gas: " << name << std::endl;
        for (uint16_t i = 0; i < stateTree.size(); ++i) {
            stateTree[i]->print();
        }
    }
};

struct EedfGas : public Gas<Boltzmann> {
    void someExtraFunction() {
        std::cout << "I am an EedfGas." << std::endl;
    }

    explicit EedfGas(const std::string &name)
            : Gas<Boltzmann>(name) {}
};

struct ChemGas : public Gas<Chemistry> {
    void someExtraFunction() {
        std::cout << "I am a ChemGas." << std::endl;
    }

    explicit ChemGas(const std::string &name)
            : Gas<Chemistry>(name) {}
};

#endif //LOKI_CPP_GAS_H
