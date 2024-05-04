#ifndef LOKI_CPP_GASMIXTUREBASE_H
#define LOKI_CPP_GASMIXTUREBASE_H

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
class GasMixture
{
  public:
    using Gases = std::vector<std::unique_ptr<Gas>>;
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
    /** Loads the data concerning a single property of the states (energy, statistical
     *  weight, or population) from a vector of entries as supplied by the setup object.
     *  First it determines whether the current entry requires loading by direct value,
     *  file or function, then it acts accordingly. The property that needs to be set
     *  is defined by the propertyType variable. Furthermore, it also needs a pointer
     *  to the WorkingConditions structure to access the argument map (which maps its
     *  member variables to names by which they are addressed in the input files).
     */
    void loadStateProperty(const std::filesystem::path &basePath, const std::vector<std::string> &entryVector, 
                           StatePropertyType propertyType, const WorkingConditions *workingConditions);
    /** Evaluates the densities of all states in the state tree. \todo Explain. All states of all gases?
     */
    void evaluateReducedDensities();
    /** Loads the state properties as specified in the input file. It also calls
     *  checkPopulations.
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
    Gas::State *findStateById(const std::string &stateId);
  private:
    Gas *addGas(const GasProperties& gasProps, const std::string& name);
    Gases m_gases;
    StateMap m_states;
};

} // namespace loki

#endif // LOKI_CPP_GASMIXTUREBASE_H
