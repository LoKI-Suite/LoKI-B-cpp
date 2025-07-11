/** \file
 *
 *  Declaration of the GasMixture class.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \date   May 2019
 */

#ifndef LOKI_CPP_GASMIXTURE_H
#define LOKI_CPP_GASMIXTURE_H

#include "LoKI-B/Exports.h"
#include "LoKI-B/Gas.h"
#include "LoKI-B/GasProperties.h"
#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/json.h"

#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <iostream>
#include <string>
#include <vector>

namespace loki
{

/** Gas mixture base class...
 */
class lokib_export GasMixture
{
  public:
    using Gases = std::vector<std::unique_ptr<Gas>>;
    GasMixture() = default;
    GasMixture(const GasMixture& other) = delete;
    GasMixture(GasMixture&& other) = delete;
    const GasMixture& operator=(const GasMixture& other) = delete;
    ~GasMixture();
    void print(std::ostream &os);
    /** Checks whether the sum of the gas fractions is equal to 1.
     */
    void checkGasFractions() const;
    void checkPopulations() const;
    // Vector of pointers to all the gases in the mixture.
    const Gases &gases() const { return m_gases; }
    /** Tries to add a new gas based on a given name and returns a pointer to it. If a
     *  gas with the same name already exists it returns a pointer to that gas.
     */
    Gas *ensureGas(const GasProperties& gasProps, const std::string &name);
    /** Tries to find a gas in the gas mixture specified by a supplied name. If the
     *  gas is present, a pointer to this gas is returned, otherwise a nullptr is
     *  returned.
     */
    const Gas *findGas(const std::string &name) const;
    /** See the non-constant overload of this member.
     */
    Gas *findGas(const std::string &name);
    /** Tries to find a state in the gas mixture specified by a supplied name. If the
     *  state is present, a pointer to this state is returned, otherwise a nullptr is
     *  returned.
     *
     *  Note that when entry describes a series of states (thus it contains a wildcard
     *  character) then a pointer to the first sibling is returned.
     */
    Gas::State *findState(const StateEntry &entry);
    /** Returns a container of matches. Empty (no match), one state (no owildcard)
     *  or a vector with one or more state pointers (entry has a wildcard).
     *  \todo reconsider all the state accessors. Which interfaces do we really need.
     */
    Gas::State::ChildContainer findStates(const StateEntry &entry);
    /** Evaluates the densities of all states in the state tree. \todo Explain. All states of all gases?
     */
    void evaluateReducedDensities();
    /** Loads the state properties energy, statisticalWeight and population
     *  as specified in subsections "energy", "statisticalWeight" and "population"
     *  of \a cnf. Each is handled by a call to private member loadStateProperty.
     */
    void loadStateProperties(const std::filesystem::path &basePath, const json_type &cnf, 
                             const WorkingConditions *workingConditions);

    using StateMap = std::map<std::string, Gas::State *>;
    /** This overload accepts only a reference to a StateEntry object. This is the main
     *  function to call when a new state is to be added to the mixture. It will add the
     *  corresponding gas and ancestor states if the do not yet exist, before adding the
     *  state itself. If the state already exists, then a pointer to this state is
     *  returned instead.
     *
     *  Note that this container only contains states that are actually present in the
     *  mixture, not their parents that are crated implicitly.
     */
    Gas::State *ensureState(const GasProperties& gasProps, const StateEntry &entry);
    /// \todo Add a constant overload of GasMixture::findStateById
    Gas::State *findStateById(const std::string &stateId);
  private:
    /** propEntry should be an object that can contain:
     *  \verbatim
           {
             "<stateID>": {
               "type": "constant",
               "value": <value>
             }
           } \endverbatim

     *  or
     *  \verbatim
           {
             "<stateID>": {
               "type": "function",
               "name": <funcname>,
               "arguments": [ <arguments> ]
             }
           } \endverbatim
     *     "<arguments>" is an array. Each argument is either a parameter name
     *     (a string) or a direct value (a double). It can be empty or be
     *     emitted entirely is the function has no arguments.
     *
     *  First it determines whether the current entry requires loading by direct value,
     *  file or function, then it acts accordingly. The property that needs to be set
     *  is defined by the propertyType variable. Furthermore, it also needs a pointer
     *  to the WorkingConditions structure to access the argument map (which maps its
     *  member variables to names by which they are addressed in the input files).
     */
    void loadStatePropertyEntry(const std::string& state_id, const json_type& propEntry,
                           StatePropertyType propertyType, const WorkingConditions *workingConditions);
    /** Sets a property (energy, statistical weight, population) of selected states
     *  as specified in json array \a stateProp. The  members of this array must be
     *  objects of the form "{ "files": [ <filenames>] }", or of the form that is
     *  expected for the propEntry argument of member \c loadStatePropertyEntry.
     *
     *  If files specification is found, a JSON object with the same structure
     *  as argument \a stateProp is created for each file, and loadStateProperty
     *  is called recursively with that JSON object as argument (as if that file
     *  were included in the place of the { "files: <filenames> } member). If the
     *  filename has extension ".json", it will be used as is. For other files
     *  the JSON object is created via a call to to readLegacyStatePropertyFile,
     *  which converts the legacy LoKI-B "*.in" file format to JSON.
     *
     *  For every "value" or "function" node, loadStatePropertyEntry is called to
     *  set the values of the property of the selected states.
     *
     *  Argument \a basePath is needed to resolve the location of included files with
     *  a relative path name. Argument \a propertyType indices which property (data
     *  member) of the states is being configures, the \a workingConditions are needed
     *  to translate parameter names (that can be used as function arguments) into their
     *  numerical values.
     *
     *  \sa loadStatePropertyEntry
     *  \sa loadStateProperties
     *  \sa statePropertyFile
     *
     *  \todo There is the risk of infinite recursion if a file refers to that
     *        same file (possibly indirectly).
     */
    void loadStateProperty(const std::filesystem::path &basePath, const json_type& stateProp,
                           StatePropertyType propertyType, const WorkingConditions *workingConditions);
    Gas *addGas(const GasProperties& gasProps, const std::string& name);
    Gases m_gases;
    StateMap m_states;
};

} // namespace loki

#endif // LOKI_CPP_GASMIXTURE_H
