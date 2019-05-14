//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_STATE_H
#define LOKI_CPP_STATE_H

#include "Enumeration.h"

#include <string>

namespace loki {
    class State {
    protected:
        Enumeration::StateType type;
        std::string e;
        uint16_t charge, v, J;
        double energy,
               statisticalWeight,
               population,
               density;

    public:
        State() = delete;
    };
}

#endif //LOKI_CPP_STATE_H
