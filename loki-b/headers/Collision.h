//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_COLLISION_H
#define LOKI_CPP_COLLISION_H

#include "Enumeration.h"
#include "EedfState.h"
#include "Grid.h"

#include <vector>

namespace loki {
    class Collision {
        Enumeration::CollisionType type;

        // isExtra might not be necessary (since all these
        // collisions will be in another array anyway).
        bool isExtra = false, isReverse = false;
        double threshold = 0.;

        State * target;
        std::vector<State *> products;
        std::vector<uint16_t> stoiCoeff;

        Grid * energyGrid;

        // TODO: decide whether to abstract RawCrossSection and
        //  CrossSection into separate classes.

        std::vector<std::pair<double, double>> rawCrossSection;
        std::vector<double> crossSection;

        // Rate Coefficient variables
    };
}


#endif //LOKI_CPP_COLLISION_H
