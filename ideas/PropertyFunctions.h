#ifndef LOKI_CPP_IDEAS_PROPERTY_FUNCTIONS_H
#define LOKI_CPP_IDEAS_PROPERTY_FUNCTIONS_H

#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "ideas/StateContainer.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Enumeration.h"

namespace loki::experimental{

/** In this namespace the user can define all their property functions. Each property function has the same
 * structure.
 *
 * The arguments are as follows.
 *
 *  1. A vector of pointers to the states (typename Trait<TraitType>::State) of which to set the property.
 *  2. A vector of doubles representing the arguments to be used by the property function. E.g. in the
 *     case of boltzmannPopulation this vector contains a single argument, the temperature.
 *  3. A StatePropertyType indicating whether the function is supposed to set either the energy, population
 *     or statistical weight of the states.
 *
 * The user is responsible to provide any checks on the number of arguments provided, the validity of the
 * passed StateProperty etc. Furthermore, when adding a function, the user should also add it to the
 * callByName function, which links the function to a string containing its name. When this is done, it is
 * possible to address the function in a setup file to set state properties.
 */
namespace PropertyFunctions {

/** Accepts a pointer to a state, a double value and a StatePropertyType. Based on the value of the
 *  StatePropertyType, it will assign the value to either the population, statisticalWeight or
 *  population member variables of the supplied state object.
 */
template <class StateT>
inline void setStateProperty(StateT& state, const double &value, loki::StatePropertyType type)
{
    switch (type) {
        case loki::StatePropertyType::energy:
            state.energy() = value;
            break;
        case loki::StatePropertyType::statisticalWeight:
            state.statisticalWeight() = value;
            break;
        case loki::StatePropertyType::population:
            state.population() = value;
            break;
    }
}

template <class StateT>
inline void constantValue(const ChildContainer<StateT> &states,
                   double value, loki::StatePropertyType type) {
//            if (arguments.size() != 1)
//                Log<NumArgumentsError>::Error("constantValue");
//
//            const double &value = arguments[0];
    for (auto& state : states)
    {
        setStateProperty(*state, value, type);
    }
}

template <class StateT>
inline void boltzmannPopulation(const ChildContainer<StateT> &states, const std::vector<double> &arguments)
{
    if (arguments.size() != 1)
        Log<NumArgumentsError>::Error("boltzmannPopulation");

    const double &temp = arguments[0];
    double norm = 0.;

    for (auto& state : states) {
        state->population() = state->statisticalWeight() * std::exp(-state->energy() / (Constant::kBeV * temp));
        norm += state->population();
    }
    for (auto& state : states) {
        state->population() /= norm;
    }
}

inline void harmonicOscillatorEnergy(const ChildContainer<StateVibrational> &states,
                              const std::vector<double> &arguments)
{
    if (!arguments.empty())
        Log<NumArgumentsError>::Error("harmonicOscillatorEnergy");

    const Gas& gas = (*states.begin())->getGas();
    if (gas.harmonicFrequency < 0)
        Log<Message>::Error("Cannot find harmonicFrequency of the gas " + gas.name() +
                            " to evaluate state energies.");

    for (auto& state : states)
    {
        double v;

        if (!Parse::getValue(state->v(), v))
            Log<Message>::Error("Non numerical vib level (" + state->v() +
                                ") when trying to assign harmonic oscillator energy.");

        state->energy() = Constant::plankReducedInEv * gas.harmonicFrequency * (v + .5);
    }
}

inline void rotationalDegeneracy_N2(const ChildContainer<StateRotational> &states,
                      const std::vector<double> &arguments)
{
    if (!arguments.empty())
        Log<NumArgumentsError>::Error("rotationalDegeneracy_N2");

    for (auto& state : states)
    {
        double J;

        if (!Parse::getValue(state->J(), J))
            Log<Message>::Error("Non numerical rot level (" + state->J() +
                                ") when trying to assign rigid rotor energy.");

        state->statisticalWeight() = 3*(1. + .5 * (1. + std::pow(-1., J))) * (2. * J + 1.);
    }
}

inline void rotationalDegeneracy(const ChildContainer<StateRotational>& states,
                             const std::vector<double> &arguments)
{
    if (!arguments.empty())
        Log<NumArgumentsError>::Error("rotationalDegeneracy");

    for (auto& state : states)
    {
        double J;

        if (!Parse::getValue(state->J(), J))
            Log<Message>::Error("Non numerical rot level (" + state->J() +
                                ") when trying to assign rigid rotor energy.");

        state->statisticalWeight() = 2 * J + 1;
    }
}


inline void rigidRotorEnergy(const ChildContainer<StateRotational>&states,
                      const std::vector<double> &arguments)
{
    if (!arguments.empty())
        Log<NumArgumentsError>::Error("rigidRotorEnergy");

    if (states.empty())
        Log<Message>::Error("rigidRotorEnergy: empty rotational state list.");
    const Gas& gas = (*states.begin())->getGas();
    if (gas.rotationalConstant < 0)
        Log<Message>::Error("Cannot find rotationalConstant of the gas " + gas.name() +
                            " to evaluate state energies.");

    for (auto& state : states)
    {
        double J;

        if (!Parse::getValue(state->J(), J))
            Log<Message>::Error("Non numerical rot level (" + state->J() +
                                ") when trying to assign rigid rotor energy.");

        state->energy() = gas.rotationalConstant * J * (J + 1.);
    }

#if 0
inline void callByName(const std::string &name, const std::vector<GasBase::StateBase*> &states,
                const std::vector<double> &arguments, StatePropertyType type) {

    if (name == "boltzmannPopulation")
        boltzmannPopulation(states, arguments, type);
    else if (name == "harmonicOscillatorEnergy")
        harmonicOscillatorEnergy(states, arguments, type);
    else if (name == "rigidRotorEnergy")
        rigidRotorEnergy(states, arguments, type);
    else if (name == "rotationalDegeneracy_N2")
        rotationalDegeneracy_N2(states, arguments, type);
    else if (name == "rotationalDegeneracy")
        rotationalDegeneracy(states, arguments, type);
    else
        Log<PropertyFunctionError>::Error(name);
}
#endif

}


} // namespace PropertyFunctions
} // namespace loki::experimental

#endif // LOKI_CPP_IDEAS_PROPERTY_FUNCTIONS_H
