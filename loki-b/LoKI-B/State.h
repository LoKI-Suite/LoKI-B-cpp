//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_STATE_H
#define LOKI_CPP_STATE_H

#include "InputStructures.h"
#include "Enumeration.h"
#include "Traits.h"
#include "GasBase.h"
#include "EedfGas.h"
#include "Log.h"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

namespace loki
{

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
template <typename TraitType>
class State : public GasBase::StateBase
{
protected:
    /** Instantiates a State object as presented by a StateEntry object. Take note that
     * the type of the state is inferred from its parent and not from the StateEntry.
     * This is done to allow one StateEntry object to be used to create a state and its
     * ancestors if they do not yet exist.
     */
    State(const StateEntry &entry,
          typename Trait<TraitType>::Gas *gas,
          typename Trait<TraitType>::State *parent)
    : StateBase(entry,parent ? static_cast<Enumeration::StateType>(parent->type + 1) : electronic, gas, parent)
    {
    }
    State(Enumeration::StateType type, typename Trait<TraitType>::Gas *gas,
          std::string e, std::string v, std::string J, std::string charge)
    : StateBase(type,e,v,J,charge,gas)
    {
    }
public:
    typename Trait<TraitType>::State* parent() { return static_cast<typename Trait<TraitType>::State*>(GasBase::StateBase::parent()); }
    typename Trait<TraitType>::Gas& gas() { return static_cast<typename Trait<TraitType>::Gas&>(GasBase::StateBase::gas()); }

    const std::vector<typename Trait<TraitType>::State *>& children() const { return m_children; }

    /* Returns a vector containing the sibling states of the current state.
     * In case this state is electronic, this vector is obtained via the
     * gas, otherwise it is obtained via the parent.
     */
    const std::vector<typename Trait<TraitType>::State *> &siblings() const
    {
        if (type == electronic)
        {
            return charge.empty() ? gas()->stateTree : gas()->ionicStates;
        }
        return parent()->m_children;
    }
    typename Trait<TraitType>::State* add_child(typename Trait<TraitType>::State* state)
    {
        StateBase::add_child(state);
        return m_children.emplace_back(state);
    }

    /* Allows to search the children in order to find a state that is equal
     * to or an ancestor of the state as described by the passed StateEntry
     * object. If this state is present, it will be returned, otherwise a
     * null pointer is returned.
     */
    typename Trait<TraitType>::State *
    find(const StateEntry &entry)
    {
        return static_cast<typename Trait<TraitType>::State*>(StateBase::find(entry));
    }

    private:

    std::vector<typename Trait<TraitType>::State *> m_children;
};

} // namespace loki

#endif //LOKI_CPP_STATE_H
