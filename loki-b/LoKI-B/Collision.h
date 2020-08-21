//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_COLLISION_H
#define LOKI_CPP_COLLISION_H

#include "Enumeration.h"

#include <vector>

namespace loki {

    /** This class is a base class to the EedfCollision and future ChemCollision (Reaction) classes.
     *  It is templated to allow the use of trait classes (classes that define types).
     *
     *  Collision<Boltzmann> stores a pointer to a single EedfState target and a vector of pointers
     *  to EedfState products.
     *  Collision<Chemistry> stores one vector of pointers to EedfStates for the LHS and one for the
     *  RHS.
     */
    template<typename StateT>
    class Collision
    {
    public:
        using State = StateT;

        virtual ~Collision() {}

        const CollisionType type;
        const bool isReverse;

        using StateVector = std::vector<State*>;
protected:
        const StateVector target;
public:

        // Vector of pointers to the product states.
        const StateVector products;

        // Vector of stoichiometric coefficents belonging to the products of the collsion.
        const std::vector<uint16_t> productStoiCoeff;

    protected:

        /** Initializes the members of the Collision class, the vectors are initialized
         *  using move semantics.
         */
        Collision(CollisionType type,
                StateVector reactants,
                StateVector products,
                std::vector<uint16_t> &stoiCoeff, bool isReverse)
        : type(type),
        isReverse(isReverse),
        target(std::move(reactants)),
        products(std::move(products)),
        productStoiCoeff(std::move(stoiCoeff))
        {
        }
    };
}


#endif //LOKI_CPP_COLLISION_H
