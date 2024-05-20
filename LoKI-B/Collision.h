/** \file
 *
 *  Interfaces of LoKI-B's code for collisions.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   2 May 2019
 */


#ifndef LOKI_CPP_COLLISION_H
#define LOKI_CPP_COLLISION_H

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Gas.h"

#include <map>
#include <vector>

namespace loki
{

/** This class is a base class to the EedfCollision and future ChemCollision (Reaction) classes.
 */
class Collision
{
  public:
    using State = Gas::State;
    using StateVector = std::vector<const State *>;
    using CoeffVector = std::vector<uint16_t>;

    virtual ~Collision()
    {
    }

    CollisionType type() const { return m_type; }
    bool isReverse() const { return m_isReverse; }

    /** Compare the stoichiometry and type of the reaction.
     */
    bool is_same_as(const Collision &other) const
    {
        return m_type == other.m_type && m_lhs == other.m_lhs && m_rhs == other.m_rhs;
    }
protected:
    /** Initializes the members of the Collision class, the vectors are initialized
     *  using move semantics.
     */
    Collision(CollisionType type, const StateVector &lhsStates, const CoeffVector &lhsCoeffs,
              const StateVector &rhsStates, const CoeffVector &rhsCoeffs, bool isReverse)
        : m_type(type), m_isReverse(isReverse), m_lhs(makeStoichMap(lhsStates, lhsCoeffs)),
          m_rhs(makeStoichMap(rhsStates, rhsCoeffs))
    {
    }
private:
    using StoichArray = std::map<const State *, uint16_t>;
    static StoichArray makeStoichMap(const StateVector &states, const CoeffVector &coeffs)
    {
        StoichArray sm;
        assert(states.size() == coeffs.size());
        for (unsigned i = 0; i != states.size(); ++i)
        {
            sm[states[i]] += coeffs[i];
        }
        return sm;
    }
    const CollisionType m_type;
    const bool m_isReverse;
    // Left-hand side-map of State pointers and stoichiometric coefficients
    const StoichArray m_lhs;
    // Right-hand side-map of State pointers and stoichiometric coefficients
    const StoichArray m_rhs;
};

} // namespace loki

#endif // LOKI_CPP_COLLISION_H
