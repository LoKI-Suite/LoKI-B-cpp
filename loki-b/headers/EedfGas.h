//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_EEDFGAS_H
#define LOKI_CPP_EEDFGAS_H

#include "Gas.h"
#include "EedfState.h"
#include "EedfCollision.h"

#include <vector>

namespace loki {
    class EedfGas : public Gas<Boltzmann> {

        // We need to store the collisions per Gas since we need to calculate
        // the mass ratio when evaluating the total and elastic cross-sections.
        std::vector<EedfCollision *> collisions, extraCollisions;

    public:
        // TODO: effectivePopulations -> what type?
        double OPBParameter = 0.;

        explicit EedfGas(const std::string &name);

        ~EedfGas();

        void addCollision(EedfCollision * collision, bool isExtra);
    };
}


#endif //LOKI_CPP_EEDFGAS_H
