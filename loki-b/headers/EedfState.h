//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFSTATE_H
#define LOKI_CPP_EEDFSTATE_H

#include "State.h"
#include "Traits.h"

#include <vector>

namespace loki {

    class EedfState : public State<Boltzmann> {
        bool isTarget{false};

        // Are these needed?
        //std::vector<Collision *> collisions, extraCollisions;

    public:
        EedfState(const StateEntry &entry, EedfGas *gas, EedfState * parent = nullptr);
        ~EedfState() override = default;

        EedfState(const EedfState &other) = delete;
    };
}


#endif //LOKI_CPP_EEDFSTATE_H
