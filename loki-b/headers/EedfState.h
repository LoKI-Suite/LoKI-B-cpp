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

        std::vector<EedfCollision *> collisions, extraCollisions;

    public:
        EedfState(const StateEntry &entry, EedfGas *gas, EedfState * parent = nullptr);
        EedfState(StateType type, EedfGas *gas, const std::string &e, const std::string &v,
                const std::string &J, const std::string &charge);

        EedfState(const EedfState &other) = delete;

        ~EedfState() override = default;

        bool hasCollision(EedfCollision * collision, bool isExtra);
        void addCollision(EedfCollision * collision, bool isExtra);

        friend std::ostream &operator<<(std::ostream &os, const EedfState &state);
    };
}


#endif //LOKI_CPP_EEDFSTATE_H
