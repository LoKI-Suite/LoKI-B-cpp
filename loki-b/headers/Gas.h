//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_GAS_H
#define LOKI_CPP_GAS_H

#include "Traits.h"

#include <iostream>
#include <string>
#include <vector>

namespace loki {
    template <typename TraitType>
    class Gas {
    protected:
        const std::string name;
        double mass{0.},
               harmonicFrequency{0.},
               anharmonicFrequency{0.},
               rotationalConstant{0.},
               electricDipoleMoment{0.},
               electricQuadrupoleMoment{0.},
               polarizability{0.},
               fraction{0.};

        std::vector<typename Trait<TraitType>::State *> states;
        std::vector<typename Trait<TraitType>::State *> stateTree;

        explicit Gas(std::string name) : name(std::move(name)) {}

    public:
        bool operator==(const std::string &otherName) {
            return name == otherName;
        }

        void print() const {
            std::cout << "Gas: " << name << std::endl;
            for (uint16_t i = 0; i < stateTree.size(); ++i) {
                stateTree[i]->print();
            }
        }

        ~Gas() {
            for (auto *state : states)
                delete state;
        }
    };
}


#endif //LOKI_CPP_GAS_H
