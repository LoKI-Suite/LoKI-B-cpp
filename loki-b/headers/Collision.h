//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_COLLISION_H
#define LOKI_CPP_COLLISION_H

#include "Enumeration.h"
#include "EedfState.h"
#include "Grid.h"
#include "CrossSection.h"
#include "Traits.h"

#include <vector>

namespace loki {
    template<typename TraitType>
    class Collision {
        Enumeration::CollisionType type;

        // isExtra might not be necessary (since all these
        // collisions will be in another array anyway).
        bool isExtra = false, isReverse = false;

        State<typename Trait<TraitType>::State> * target;
        std::vector<typename Trait<TraitType>::State *> products;
        std::vector<uint16_t> stoiCoeff;

        Grid * energyGrid;

        // The raw cross section data an threshold is stored in
        // the CrossSection object

        CrossSection crossSection;

        // Rate Coefficient variables should be here
    };
}


#endif //LOKI_CPP_COLLISION_H
