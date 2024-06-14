/** \file
 *
 *  Declaration of a class that manages characteristic numbers of a state
 *  and code to retrieve these numbers from a character string.
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
 *  \author Daan Boer and Jan van Dijk
 *  \date   23 May 2019
 */

#ifndef LOKI_CPP_STATEENTRY_H
#define LOKI_CPP_STATEENTRY_H

#include <filesystem>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/json.h"
#include "LoKI-B/Exports.h"

namespace loki
{

// using namespace Enumeration;

class lokib_export StateEntry
{
  public:
    StateEntry();
    StateEntry(const std::string &id, StateType level, const std::string &gasName, const std::string &charge,
               const std::string &e, const std::string &v, const std::string &J);
    static StateEntry electronEntry();
    bool hasWildCard() const;

    /** \todo Rethink m_id, or at least document this. id should only be associated with
     *        the final state that is created, not with intermediate parents. This is not
     *        a problem now, because id is only the in the map that is managed by the mixture
     *        in which only this final result is stored, but still: can an id be assigned
     *        also to intermediate states in a reliable way?
     */
    const std::string m_id;
    const StateType m_level;
    const std::string m_gasName, m_charge, m_e, m_v, m_J;
};

lokib_export std::ostream& operator<<(std::ostream &os, const StateEntry &entry);

/* Accepts a string containing the LHS or RHS of a collision equation. The states
 * in this expression are then parsed into a vector of StateEntry objects, which
 * is passed by reference. Furthermore, the user can supply a pointer to a vector
 * in which the stoichiometric coefficients of the states in this collision are
 * then stored. If a null pointer is passed, these coefficients are not stored.
 */
lokib_export void entriesFromString(const std::string stateString, std::vector<StateEntry> &entries,
                       std::vector<uint16_t> *stoiCoeff);

/// \todo Make this a StateEntry constructor
lokib_export StateEntry entryFromJSON(const std::string& id, const json_type &cnf);

/** Extracts a StateEntry object from a given string and returns it. Note that this function
 *  is specifically used when loading state properties, since then the states can contain
 *  wild card characters.
 */
lokib_export StateEntry propertyStateFromString(const std::string &propertyString);

} // namespace loki

#endif // LOKI_CPP_STATEENTRY_H
