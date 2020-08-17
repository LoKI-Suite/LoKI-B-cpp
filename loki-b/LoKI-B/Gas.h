//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_GAS_H
#define LOKI_CPP_GAS_H

#include "GasBase.h"
#include "Log.h"

#include <string>
#include <vector>
#include <type_traits>

namespace loki {
    /** This class acts as a base class for both the EedfGas and the future ChemGas. It
     *  is templated to allow the use of trait classes (classes that define types). The
     *  constructor of this class is protected such that the user cannot instantiate it
     *  as is. However, other classes can inherit from this class. In this case we can
     *  either inherit Gas<Boltzmann> or Gas<Chemistry>. The difference here is that
     *  the Gas<Boltzmann> class has vectors to pointers of EedfStates whereas
     *  Gas<Chemistry> has vectors to pointers of ChemStates.
     *
     *  The reason that this class is introduced in this format is that an EedfGas and
     *  a ChemGas share a lot of properties. On top of that they share the same structure
     *  i.e. they both contain vectors to states, however the type of the states differs
     *  between the two, which is fixed by introducing the trait classes. In this way
     *  we can implement some of the functionality, such as adding states to a gas,
     *  only once instead of separately for both classes.
     */

    template <typename GasType>
    class Gas : public GasBase {
    public:
        class State;

        State& get_root() { return static_cast<State&>(*GasBase::m_root); }
        const State& get_root() const { return static_cast<State&>(*GasBase::m_root); }

        /// \todo Needed?
        State* findState(const StateEntry& entry) {
            return static_cast<State*>(GasBase::findState(entry));
        }
        State* ensureState(const StateEntry &entry)
        {
            auto *state = get_root().ensure_state(entry);
            // state is now a 'charge state' (container)
            for (uint8_t lvl = charge; lvl < entry.level; ++lvl) {
                state = state->ensure_state(entry);
            }
            return state;
        }
    protected:
        explicit Gas(std::string name)
        : GasBase(name)
        {
            static_assert(std::is_base_of_v<Gas,GasType>);
            this->m_root.reset(new State(static_cast<GasType*>(this)));
        }
    };

/** State is a base class to the EedfState and ChemState classes. It is templated
 * to allow for the use of trait classes (classes that define types). The
 * constructor of this class is protected such that the user cannot instantiate
 * this class as is. However, other classes can be derived from this class. They
 * can either inherit from State<Boltzmann> or State<Chemistry>. The difference
 * here is that a State<Boltzmann> has EedfStates as its parent and children,
 * and a pointer to and EedfGas, whereas a State<Chemistry> has ChemStates as
 * its parent and children and a pointer to a ChemGas.
 *
 * The idea behind the implementation of the State base class is similar to that
 * of the Gas base class.
 */
template <typename GasType>
class Gas<GasType>::State : public GasBase::StateBase
{
public:
    using Gas = GasType;
    /** Instantiates a State object as presented by a StateEntry object. Take note that
     * the type of the state is inferred from its parent and not from the StateEntry.
     * This is done to allow one StateEntry object to be used to create a state and its
     * ancestors if they do not yet exist.
     */
    State(const StateEntry &entry,
          Gas *gas,
          State &parent)
    : StateBase(entry, gas, parent)
    {
    }
    State(Gas *gas)
    : StateBase(gas)
    {}

    State* parent() { return static_cast<State*>(GasBase::StateBase::parent()); }
    Gas& gas() { return static_cast<Gas&>(GasBase::StateBase::gas()); }

    const std::vector<State *>& children() const { return m_children; }

    /* Returns a vector containing the sibling states of the current state.
     * In case this state is electronic, this vector is obtained via the
     * gas, otherwise it is obtained via the parent.
     */
    const std::vector<State *> &siblings() const
    {
        if (!parent())
        {
            throw std::runtime_error("Cannot call siblings on the root node.");
        }
        return parent()->m_children;
    }
    State* ensure_state(const StateEntry &entry)
    {
        State* state = find(entry);
        if (state)
        {
            return state;
        }
        state = new State(entry, &gas(), *this);
        StateBase::add_child(state);
        m_children.emplace_back(state);
        return state;
    }

    /* Allows to search the children in order to find a state that is equal
     * to or an ancestor of the state as described by the passed StateEntry
     * object. If this state is present, it will be returned, otherwise a
     * null pointer is returned.
     */
    State* find(const StateEntry &entry)
    {
        return static_cast<State*>(StateBase::find(entry));
    }

    private:

    std::vector<State *> m_children;
};

} // namespace loki

#endif //LOKI_CPP_GAS_H
