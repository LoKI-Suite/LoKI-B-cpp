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
    /* -- Collision --
     * This class is a base class to the EedfCollision and future ChemCollision (Reaction) classes.
     * It is templated to allow the use of trait classes (classes that define types).
     *
     * Collision<Boltzmann> stores a pointer to a single EedfState target and a vector of pointers
     * to EedfState products.
     * Collision<Chemistry> stores one vector of pointers to EedfStates for the LHS and one for the
     * RHS.
     */

    template<typename TraitType>
    class Collision {
    public:
        Enumeration::CollisionType type;
        bool isReverse = false;

        // In the case that TraitType = Boltzmann, this will be a pointer to an EedfState.
        // When TraitType = Chemistry, this is a vector of pointers to Chemstates.
        typename Trait<TraitType>::Reactants target;

        // Vector of pointers to the product states.
        std::vector<typename Trait<TraitType>::State *> products;

        // Vector of stoichiometric coefficents belonging to the products of the collsion.
        std::vector<uint16_t> productStoiCoeff;

    protected:

        /* -- Constructor --
         * Initializes the members of the Collision class, the vectors are initialized
         * using move semantics.
         */
        Collision(Enumeration::CollisionType type,
                typename Trait<TraitType>::Reactants reactants,
                std::vector<typename Trait<TraitType>::State *> products,
                std::vector<uint16_t> &stoiCoeff, bool isReverse)
                : type(type), target(std::move(reactants)), products(std::move(products)),
                  productStoiCoeff(std::move(stoiCoeff)), isReverse(isReverse) {}
    };
}


#endif //LOKI_CPP_COLLISION_H
