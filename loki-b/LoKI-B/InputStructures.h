//
// Created by daan on 23-5-19.
//

#ifndef LOKI_CPP_INPUTSTRUCTURES_H
#define LOKI_CPP_INPUTSTRUCTURES_H

#include <ostream>
#include <cstdint>
#include <string>
#include <vector>

#include "Enumeration.h"
#include "LoKI-B/json.h"

namespace loki {
    using namespace Enumeration;


    struct StateEntry {
        StateType level;
        std::string charge, gasName, e, v, J;

        StateEntry() : level(none) {}

        StateEntry(StateType level, const std::string &gasName, const std::string &charge,
                   const std::string &e, const std::string &v, const std::string &J)
                : level(level), charge(charge), gasName(gasName),
                  e(e), v(v), J(J) {}

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

        friend std::ostream &operator<<(std::ostream &os, const StateEntry &entry) {
            os << entry.gasName << '(';

            if (!entry.charge.empty()) os << entry.charge << ',';

            os << entry.e << ",v=" << entry.v << ",J=" << entry.J << ')';

            return os;
        }
    };

}

#endif //LOKI_CPP_INPUTSTRUCTURES_H
