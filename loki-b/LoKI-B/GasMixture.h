//
// Created by daan on 14-5-19.
//

#ifndef LOKI_CPP_GASMIXTURE_H
#define LOKI_CPP_GASMIXTURE_H

// DONE: Fill the Gas, State and Collision structures.
// DONE: Implement the templated approach to gases and states
// DONE: Think of how we can parse the LXCat data. Especially for states, do we need an
//  auxiliary structure like StateEntry or can we immediately create the state?
//  If not then alter the GasMixture class (replace StateEntry by State).
// DONE: Implement LXCat file parsing.
// DONE: Implement Property file / function / value parsing.

#include "LoKI-B/GasMixtureBase.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/json.h"

#include <fstream>
#include <map>
#include <vector>

namespace loki
{

/** The GasMixture class acts as a base class to the EedfGasMixture and future
 *  ChemGasMixture classes. This class is templated to allow the use of trait
 *  classes (classes that define types). Its constructor is protected such
 *  that the class cannot be instantiated as is, although other classes can be
 *  derived from this class. They can either inherit from GasMixture<Boltzmann>
 *  or GasMixture<Chemistry>. The idea is that the first handles EedfGases,
 *  States and Collisions, whereas the latter handles their Chemistry
 *  equivalents.
 */
template <typename GasT>
class GasMixture : public GasMixtureBase
{
  public:
    using Gas = GasT;
    using State = typename Gas::State;
    // Vector of pointers to all the gases in the mixture.
    const std::vector<Gas *> &gases() const
    {
        return m_gases;
    }
    /// \todo Still needed?
    Gas *findGas(const std::string &name)
    {
        return static_cast<Gas *>(GasMixtureBase::findGas(name));
    }
    /// \todo Still needed?
    State *findState(const StateEntry &entry)
    {
        return static_cast<State *>(GasMixtureBase::findState(entry));
    }
    GasBase::StateBase *findStateById(const std::string &stateId)
    {
        typename StateMap::iterator it = m_states.find(stateId);
        return it == m_states.end() ? nullptr : it->second;
    }
    /** Tries to add a new gas based on a given name and returns a pointer to it. If a
     *  gas with the same name already exists it returns a pointer to that gas.
     */
    Gas *ensureGas(const std::string &name)
    {
        Gas *gas = findGas(name);
        if (!gas)
        {
            gas = new Gas(name);
            GasMixtureBase::addGas(gas);
            m_gases.emplace_back(gas);
        }
        return gas;
    }

    using StateMap = std::map<std::string, State *>;
    /** This overload accepts only a reference to a StateEntry object. This is the main
     *  function to call when a new state is to be added to the mixture. It will add the
     *  corresponding gas and ancestor states if the do not yet exist, before adding the
     *  state itself. If the state already exists, then a pointer to this state is
     *  returned instead.
     *
     *  Note tha this container only contains states that are actually present in the
     *  mixture, not their parents that are crated implicitly.
     */
    State *ensureState(const StateEntry &entry)
    {
        typename StateMap::iterator it = m_states.find(entry.m_id);
        if (it != m_states.end())
        {
            return it->second;
        }
        State *state = ensureGas(entry.gasName)->ensureState(entry);
        m_states[entry.m_id] = state;
        return state;
    }

  private:
    // Vector of pointers to all the gases in the mixture.
    std::vector<Gas *> m_gases;
    StateMap m_states;
};

} // namespace loki

#endif // LOKI_CPP_GASMIXTURE_H
