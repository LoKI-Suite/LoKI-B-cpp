//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFGAS_H
#define LOKI_CPP_EEDFGAS_H

#include "Gas.h"
#include "EedfState.h"

#include <vector>

namespace loki {
    class EedfGas : public Gas {
        double OPBParameter;
        std::vector<EedfState> states;

        // Is this necessary? Since we already have a big array of
        // std::vector<Collision *> collisions, extraCollisions;

        // TODO: effectivePopulations -> what type?
    };
}


#endif //LOKI_CPP_EEDFGAS_H
