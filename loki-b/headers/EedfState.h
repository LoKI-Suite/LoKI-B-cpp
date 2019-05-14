//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFSTATE_H
#define LOKI_CPP_EEDFSTATE_H

#include "State.h"

#include <vector>

namespace loki {
    // Forward declaration of EedfGas
    class EedfGas;

    class EedfState : public State {
        bool isTarget{false};

        EedfGas * gas;
        EedfState * parent;

        std::vector<EedfState *> children;

        // Are these needed?
        //std::vector<Collision *> collisions, extraCollisions;

    public:
        EedfState() = default;
        ~EedfState() = default;

        EedfState(const EedfState &other) = delete;

        bool operator==(const EedfState &other);
    };
}


#endif //LOKI_CPP_EEDFSTATE_H
