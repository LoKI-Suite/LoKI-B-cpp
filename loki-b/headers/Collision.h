//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_COLLISION_H
#define LOKI_CPP_COLLISION_H

#include "Enumeration.h"
#include "EedfState.h"
#include "Traits.h"

#include <vector>

namespace loki {
    template<typename TraitType>
    class Collision {
    protected:
        Enumeration::CollisionType type;

        // isExtra might not be necessary (since all extra
        // collisions will be in another array anyway).
        bool isExtra = false, isReverse = false;

        typename Trait<TraitType>::Reactants target;
        std::vector<typename Trait<TraitType>::State *> products;
        std::vector<uint16_t> stoiCoeff;

        Collision(Enumeration::CollisionType type,
                typename Trait<TraitType>::Reactants reactants,
                std::vector<typename Trait<TraitType>::State *> products,
                std::vector<uint16_t> &stoiCoeff, bool isReverse)
                : type(type), target(std::move(reactants)), products(std::move(products)),
                  stoiCoeff(std::move(stoiCoeff)), isReverse(isReverse) {}
    };
}


#endif //LOKI_CPP_COLLISION_H
