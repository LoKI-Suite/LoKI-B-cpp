//
// Created by daan on 23-5-19.
//

#ifndef LOKI_CPP_STATEENTRY_H
#define LOKI_CPP_STATEENTRY_H

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

/** Parses state property file \a fileName into \a entries, a vector of
 *  StateEntry/double pairs.
 */
lokib_export void statePropertyFile(const std::string &fileName, std::vector<std::pair<StateEntry, double>> &entries);

} // namespace loki

#endif // LOKI_CPP_STATEENTRY_H
