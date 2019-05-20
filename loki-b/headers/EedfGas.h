//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFGAS_H
#define LOKI_CPP_EEDFGAS_H

#include "Gas.h"
#include "EedfState.h"
#include "Collision.h"

#include <vector>

namespace loki {
    class EedfGas : public Gas<Boltzmann> {
        double OPBParameter = 0.;
        std::vector<EedfState *> states;

        // We need to store the collisions per Gas since we need to calculate
        // the mass ratio when evaluating the total and elastic cross-sections.
        std::vector<Collision<Boltzmann>> collisions, extraCollisions;

        // TODO: effectivePopulations -> what type?

    public:
        explicit EedfGas(const std::string &name);
    };
}


#endif //LOKI_CPP_EEDFGAS_H
