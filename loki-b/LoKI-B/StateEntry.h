//
// Created by daan on 23-5-19.
//

#ifndef LOKI_CPP_STATEENTRY_H
#define LOKI_CPP_STATEENTRY_H

#include <ostream>
#include <string>
#include <vector>
#include <set>

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/json.h"

namespace loki {

//using namespace Enumeration;

class StateEntry
{
public:
    StateEntry();
    StateEntry(const std::string& id, StateType level, const std::string &gasName, const std::string &charge,
               const std::string &e, const std::string &v, const std::string &J);
    static StateEntry electronEntry();
    bool hasWildCard();

    const std::string m_id;
    StateType level;
    std::string charge, gasName, e, v, J;
};

std::ostream &operator<<(std::ostream &os, const StateEntry &entry);

/* Accepts a string containing the LHS or RHS of a collision equation. The states
 * in this expression are then parsed into a vector of StateEntry objects, which
 * is passed by reference. Furthermore, the user can supply a pointer to a vector
 * in which the stoichiometric coefficients of the states in this collision are
 * then stored. If a null pointer is passed, these coefficients are not stored.
 */
void entriesFromString(const std::string stateString, std::vector<StateEntry>& entries, std::vector<uint16_t>* stoiCoeff);

/// \todo Make this a StateEntry constructor
StateEntry entryFromJSON(const json_type& cnf);

/** Extracts a StateEntry object from a given string and returns it. Note that this function
 *  is specifically used when loading state properties, since then the states can contain
 *  wild card characters.
 */
StateEntry propertyStateFromString(const std::string &propertyString);

/** Parses a state property file into a vector of StateEntry, double pairs. This vector
 *  is passed by reference.
 */
bool statePropertyFile(const std::string &fileName, std::vector<std::pair<StateEntry, double>> &entries);

} // namespace loki

#endif //LOKI_CPP_STATEENTRY_H
