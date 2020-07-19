//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_STATE_H
#define LOKI_CPP_STATE_H

#include "InputStructures.h"
#include "Enumeration.h"
#include "Traits.h"
#include "Log.h"
#include "EedfGas.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

namespace loki
{
using namespace Enumeration;

/* -- State --
     * State is a base class to the EedfState and ChemState classes. It is templated
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
class State
{
protected:
    State(const StateEntry &entry,
          typename Trait<TraitType>::Gas *gas,
          typename Trait<TraitType>::State *parent);
    State(StateType type, typename Trait<TraitType>::Gas *gas,
          std::string e, std::string v, std::string J, std::string charge);

public:
    StateType type;

    typename Trait<TraitType>::State *parent{nullptr};
    typename Trait<TraitType>::Gas *gas{nullptr};

    std::string e, v, J, charge;

    double energy{-1.},
        statisticalWeight{-1.},
        population{0.},
        density{0.};

    std::vector<typename Trait<TraitType>::State *> children;

    /* -- siblings --
         * Returns a vector containing the sibling states of the current state.
         * In case this state is electronic, this vector is obtained via the
         * gas, otherwise it is obtained via the parent.
         */

    std::vector<typename Trait<TraitType>::State *> &siblings();

    /* -- >= --
         * This operator overload checks whether the current state is equal to
         * or a parent of a state as described by a given StateEntry object.
         */

    bool operator>=(const StateEntry &entry);

    template <typename T>
    friend std::ostream &operator<<(std::ostream &os, const State<T> &state);

    /* -- find --
         * Allows to search the children in order to find a state that is equal
         * to or an ancestor of the state as described by the passed StateEntry
         * object. If this state is present, it will be returned, otherwise a
         * null pointer is returned.
         */

    typename Trait<TraitType>::State *
    find(const StateEntry &entry)
    {

        auto it = std::find_if(children.begin(), children.end(),
                               [&entry](typename Trait<TraitType>::State *child) {
                                   return *child >= entry;
                               });

        if (it == children.end())
        {
            return nullptr;
        }

        return *it;
    }

    /* -- checkPopulations --
         * Verifies that the populations of all child states add up to 1. It also
         * calls the checkPopulation function on all these states to recursively check
         * populations.
         */

    void checkPopulations()
    {
        if (children.empty())
            return;

        double totalPopulation = 0.;

        for (auto *state : children)
        {
            totalPopulation += state->population;
            state->checkPopulations();
        }

        if (population == 0)
        {
            if (totalPopulation != 0)
                Log<ChildrenPopulationError>::Error(*this);
        }
        else
        {
            if (std::abs(totalPopulation - 1.) > 10. * std::numeric_limits<double>::epsilon())
                Log<ChildrenPopulationError>::Error(*this);
        }
    }

    /* -- evaluateDensity --
         * Calculates the density based on the gas fraction, the states population and populations
         * of any parent states. Then it calls evaluateDensity on all child states to recursively
         * evaluate all densities in the state tree.
         */

    void evaluateDensity();

    void printChildren() const;

    virtual ~State() = default;
};

template <typename TraitType>
State<TraitType>::State(StateType type, typename Trait<TraitType>::Gas *gas,
                        std::string e, std::string v, std::string J, std::string charge)
    : type(type), gas(gas), e(e), v(v), J(J), charge(charge) {
#if 0
    std::cout << "State::State from params: " << std::endl
	<< " type = " << static_cast<int>(type) << std::endl
	<< " gas  = " << gas << std::endl
	<< " e    = " << e << std::endl
	<< " v    = " << v << std::endl
	<< " J    = " << J << std::endl
	<< " q    = " << charge << std::endl
        << std::endl;
#endif
}

/* -- Constructor --
     * Instantiates a State object as presented by a StateEntry object. Take note that
     * the type of the state is inferred from its parent and not from the StateEntry.
     * This is done to allow one StateEntry object to be used to create a state and its
     * ancestors if they do not yet exist.
     */

template <typename TraitType>
State<TraitType>::State(const StateEntry &entry,
                        typename Trait<TraitType>::Gas *gas, typename Trait<TraitType>::State *parent)
    : parent(parent), gas(gas), e(entry.e), v(entry.v), J(entry.J), charge(entry.charge)
{
    if (parent == nullptr)
        type = electronic;
    else
        type = static_cast<StateType>(parent->type + 1);

#if 0
    std::cout << "State::State from entry: " << std::endl
	<< " type = " << static_cast<int>(type) << std::endl
	<< " gas  = " << gas << std::endl
	<< " e    = " << e << std::endl
	<< " v    = " << v << std::endl
	<< " J    = " << J << std::endl
	<< " q    = " << charge << std::endl
        << std::endl;
#endif
}

template <typename TraitType>
void State<TraitType>::printChildren() const
{
    std::string space = "  ";

    for (uint8_t i = electronic; i < type; ++i)
        space.append("  ");

    for (auto state : children)
    {
        if (state != nullptr)
        {
            std::cout << space << *state << std::endl;
            state->printChildren();
        }
    }
}

template <typename TraitType>
bool State<TraitType>::operator>=(const StateEntry &entry)
{
    {
        if (type > entry.level || charge != entry.charge)
            return false;

        if (e == entry.e)
        {
            if (type == electronic)
                return true;
        }
        else
        {
            return false;
        }

        if (v == entry.v)
        {
            if (type == vibrational)
                return true;
        }
        else
        {
            return false;
        }

        if (J == entry.J)
        {
            if (type == rotational)
                return true;
        }
        else
        {
            return false;
        }

        return true;
    }
}

template <typename TraitType>
std::vector<typename Trait<TraitType>::State *> &
State<TraitType>::siblings()
{
    if (type == electronic)
    {
        return charge.empty() ? gas->stateTree : gas->ionicStates;
    }

    return parent->children;
}

template <typename TraitType>
void State<TraitType>::evaluateDensity()
{
    if (type == electronic)
        density = population * gas->fraction;
    else
        density = population * parent->density;

    for (auto *state : children)
        state->evaluateDensity();
}

template <typename TraitType>
std::ostream &operator<<(std::ostream &os, const State<TraitType> &state)
{
    os << state.gas->name << '(';

    if (!state.charge.empty())
        os << state.charge << ',';

    os << state.e;

    if (state.type >= vibrational)
        os << ",v=" << state.v;

    if (state.type == rotational)
        os << ",J=" << state.J;

    os << ')';

    //        os << "\tpopulation: " << state.population << "\tdensity: " << state.density;

    return os;
}
} // namespace loki

#endif //LOKI_CPP_STATE_H
