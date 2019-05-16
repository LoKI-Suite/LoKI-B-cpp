//
// Created by daan on 16-5-19.
//

#ifndef LOKI_CPP_BASIC_H
#define LOKI_CPP_BASIC_H

#include <cstdint>
#include <string>

enum StateLevel : uint8_t {
    ionic,
    electronic,
    vibrational,
    rotational,
    nana
};

struct StateEntry {
    StateLevel level;
    std::string gasName, e;
    uint16_t v, J;
    int16_t charge;
};

#endif //LOKI_CPP_BASIC_H
