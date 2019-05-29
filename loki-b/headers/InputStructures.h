//
// Created by daan on 23-5-19.
//

#ifndef LOKI_CPP_INPUTSTRUCTURES_H
#define LOKI_CPP_INPUTSTRUCTURES_H

#include <cstdint>
#include <string>
#include <vector>

#include "Enumeration.h"

namespace loki {
    using namespace Enumeration;

    struct RawLXCatEntry {
        std::string reactants,
                isReverse,
                products,
                type,
                threshold;
    };

    struct StateEntry {
        StateType level;
        std::string charge, gasName, e, v, J;

        StateEntry() : level(none) {}

        StateEntry(StateType level, std::string &&gasName, std::string &&charge,
                   std::string &&e, std::string &&v, std::string &&J)
                : level(level), charge(std::move(charge)), gasName(std::move(gasName)),
                  e(std::move(e)), v(std::move(v)), J(std::move(J)) {}

        bool hasWildCard() {
            switch (level) {
                case electronic:
                    return (e == "*");
                case vibrational:
                    return (v == "*");
                case rotational:
                    return (J == "*");
                case none:
                    return false;
            }

            return false;
        }
    };

    struct CollisionEntry {
        std::vector <StateEntry> reactants, products;
        std::vector <uint16_t> stoiCoeff;
        Enumeration::CollisionType type;
        double threshold;
        bool isReverse;
//        std::vector<std::pair<double, double>> rawCrossSection;
    };
}

#endif //LOKI_CPP_INPUTSTRUCTURES_H
