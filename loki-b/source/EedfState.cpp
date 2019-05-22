//
// Created by daan on 2-5-19.
//

#include "EedfState.h"
#include "EedfGas.h"

namespace loki {


    EedfState::EedfState(const StateEntry &entry, EedfGas *gas, EedfState *parent)
        : State(entry, gas, parent) {}

    EedfState::EedfState(StateType type, EedfGas *gas, const std::string &e, const std::string &v,
                         const std::string &J, const std::string &charge)
                         : State(type, gas, e, v, J, charge) {}

    bool EedfState::hasCollision(EedfCollision *collision, bool isExtra) {
        auto &currentCollisions = (isExtra ? this->extraCollisions : this->collisions);

        return (std::find_if(currentCollisions.begin(), currentCollisions.end(),
                [collision](EedfCollision *itCollision){

                    return *itCollision == *collision;
                }) != currentCollisions.end());
    }

    void EedfState::addCollision(EedfCollision * collision, bool isExtra) {
        (isExtra ? this->extraCollisions : this->collisions).emplace_back(collision);

        // The state should only be seen as a target when it appears in a non-extra collision.
        if (!isExtra) isTarget = true;
    }

    std::ostream &operator<<(std::ostream &os, const EedfState &state) {
        os << state.gas->name << '(';

        if (!state.charge.empty()) os << state.charge << ',';

        os << state.e;

        if (state.type >= vibrational) os << ",v=" << state.v;

        if (state.type == rotational) os << ",J=" << state.J;

        os << ')';

        return os;
    }
}