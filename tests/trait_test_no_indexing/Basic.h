//
// Created by daan on 16-5-19.
//

#ifndef LOKI_CPP_BASIC_H
#define LOKI_CPP_BASIC_H

#include <cstdint>
#include <string>

// TODO: remove the 'ionic' level and store all the ionic states in a separate tree.

enum StateLevel : uint8_t {
    ionic,
    electronic,
    vibrational,
    rotational,
    none
};

struct StateEntry {
    StateLevel level;
    std::string gasName, e;
    uint16_t v, J;
    int16_t charge;
};

#endif //LOKI_CPP_BASIC_H
