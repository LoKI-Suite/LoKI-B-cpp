//
// Created by daan on 2-5-19.
//

#include "EedfState.h"
#include "EedfGas.h"

namespace loki {


    EedfState::EedfState(const StateEntry &entry, EedfGas *gas, EedfState *parent)
        : State(entry, gas, parent) {}
}