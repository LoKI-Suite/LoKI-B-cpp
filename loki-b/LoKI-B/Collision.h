//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_COLLISION_H
#define LOKI_CPP_COLLISION_H

#include "Enumeration.h"
#include "GasBase.h"

#include <vector>

namespace loki {

    /** This class is a base class to the EedfCollision and future ChemCollision (Reaction) classes.
     */
    class Collision
    {
    public:
        using State = GasBase::State;

        virtual ~Collision() {}

        const CollisionType type;
        const bool isReverse;

        using StateVector = std::vector<State*>;
        using CoeffVector = std::vector<uint16_t>;

    protected:

        /** Initializes the members of the Collision class, the vectors are initialized
         *  using move semantics.
         */
        Collision(CollisionType type,
                const StateVector& lhsStates,
                const CoeffVector& lhsCoeffs,
                const StateVector& rhsStates,
                const CoeffVector& rhsCoeffs,
                bool isReverse)
        : type(type),
        isReverse(isReverse)
#if 0
,
        m_lhsStates(lhsStates),
        m_lhsCoeffs(lhsCoeffs),
        m_rhsStates(rhsStates),
        m_rhsCoeffs(rhsCoeffs)
#endif
        {
        }
private:
#if 0
        // Vector of pointers to the reactant states.
        const StateVector m_lhsStates;
        // Vector of stoichiometric coefficents belonging to the reactants of the collsion.
        const CoeffVector m_lhsCoeffs;

        // Vector of pointers to the product states.
        const StateVector m_rhsStates;
        // Vector of stoichiometric coefficents belonging to the products of the collsion.
        const CoeffVector m_rhsCoeffs;
#endif
    };
}


#endif //LOKI_CPP_COLLISION_H
