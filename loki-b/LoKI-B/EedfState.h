//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFSTATE_H
#define LOKI_CPP_EEDFSTATE_H

#include "LoKI-B/State.h"
#include "LoKI-B/Traits.h"

#include <vector>

namespace loki {

    class EedfState : public State<Boltzmann> {
    public:
        EedfState(const StateEntry &entry, EedfGas *gas, EedfState * parent = nullptr)
            : State(entry, gas, parent)
        {}
        EedfState(StateType type, EedfGas *gas, const std::string &e, const std::string &v,
                const std::string &J, const std::string &charge)
            : State(type, gas, e, v, J, charge)
        {}

        EedfState(const EedfState &other) = delete;

        ~EedfState() override = default;
    };
}


#endif //LOKI_CPP_EEDFSTATE_H
