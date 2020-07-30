//
// Created by daan on 28-5-19.
//

#ifndef LOKI_CPP_PROPERTYFUNCTIONS_H
#define LOKI_CPP_PROPERTYFUNCTIONS_H

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Traits.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/GasBase.h"

#include <cmath>

/* -- PropertyFunctions --
 * In this namespace the user can define all their property functions. Each property function has the same
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

namespace loki::PropertyFunctions {
    using namespace Constant;

    /* -- setStateProperty --
     * Accepts a pointer to a state, a double value and a StatePropertyType. Based on the value of the
     * StatePropertyType, it will assign the value to either the population, statisticalWeight or
     * population member variables of the supplied state object.
     */

    inline void setStateProperty(GasBase::StateBase *state, const double &value, StatePropertyType type) {
        switch (type) {
            case StatePropertyType::energy:
                state->energy = value;
                break;
            case StatePropertyType::statisticalWeight:
                state->statisticalWeight = value;
                break;
            case StatePropertyType::population:
                state->population = value;
                break;
        }
    }

    inline void boltzmannPopulation(const std::vector<GasBase::StateBase*> &states,
                             const std::vector<double> &arguments, StatePropertyType type) {
        if (type != StatePropertyType::population)
            Log<WrongPropertyError>::Error("boltzmannPopulation");

        if (arguments.size() != 1)
            Log<NumArgumentsError>::Error("boltzmannPopulation");

        const double &temp = arguments[0];
        double norm = 0.;

        for (auto *state : states) {
            state->population = state->statisticalWeight * std::exp(-state->energy / (kBeV * temp));
            norm += state->population;
        }
        for (auto *state : states) {
            state->population /= norm;
        }
    }

    inline void harmonicOscillatorEnergy(const std::vector<GasBase::StateBase*> &states,
                                  const std::vector<double> &arguments, StatePropertyType type) {
        if (type != StatePropertyType::energy)
            Log<WrongPropertyError>::Error("harmonicOscillatorEnergy");

        if (!arguments.empty())
            Log<NumArgumentsError>::Error("harmonicOscillatorEnergy");

        if (states.at(0)->type != vibrational)
            Log<Message>::Error("Trying to assign harmonic oscillator energy to non-vibrational state.");

        if (states.at(0)->gas_base().harmonicFrequency < 0)
            Log<Message>::Error("Cannot find harmonicFrequency of the gas " + states.at(0)->gas_base().name +
                                " to evaluate state energies.");

        for (auto *state : states) {
            double vibLevel;

            if (!Parse::getValue(state->v, vibLevel))
                Log<Message>::Error("Non numerical vib level (" + state->v +
                                    ") when trying to assign harmonic oscillator energy.");

            state->energy = plankReducedInEv * state->gas_base().harmonicFrequency * (vibLevel + .5);
        }
    }

    inline void rigidRotorEnergy(const std::vector<GasBase::StateBase*> &states,
                          const std::vector<double> &arguments, StatePropertyType type) {
        if (type != StatePropertyType::energy)
            Log<WrongPropertyError>::Error("rigidRotorEnergy");

        if (!arguments.empty())
            Log<NumArgumentsError>::Error("rigidRotorEnergy");

        if (states.at(0)->type != rotational)
            Log<Message>::Error("Trying to assign rigid rotor energy to non-rotational state.");

        if (states.at(0)->gas_base().rotationalConstant < 0)
            Log<Message>::Error("Cannot find rotationalConstant of the gas " + states.at(0)->gas_base().name +
                                " to evaluate state energies.");

        for (auto *state : states) {
            double rotLevel;

            if (!Parse::getValue(state->J, rotLevel))
                Log<Message>::Error("Non numerical rot level (" + state->J +
                                    ") when trying to assign rigid rotor energy.");

            state->energy = state->gas_base().rotationalConstant * rotLevel * (rotLevel + 1.);
        }
    }

    inline void rotationalDegeneracy_N2(const std::vector<GasBase::StateBase*> &states,
                          const std::vector<double> &arguments, StatePropertyType type) {
        if (type != StatePropertyType::statisticalWeight)
            Log<WrongPropertyError>::Error("rotationalDegeneracy_N2");

        if (!arguments.empty())
            Log<NumArgumentsError>::Error("rotationalDegeneracy_N2");

        if (states.at(0)->type != rotational)
            Log<Message>::Error("Trying to assign rotational degeneracy to non-rotational state.");

        for (auto *state : states) {
            double J;

            if (!Parse::getValue(state->J, J))
                Log<Message>::Error("Non numerical rot level (" + state->J +
                                    ") when trying to assign rigid rotor energy.");

            state->statisticalWeight = 3*(1. + .5 * (1. + std::pow(-1., J))) * (2. * J + 1.);
        }
    }

    inline void rotationalDegeneracy(const std::vector<GasBase::StateBase*> &states,
                                 const std::vector<double> &arguments, StatePropertyType type) {
        if (type != StatePropertyType::statisticalWeight)
            Log<WrongPropertyError>::Error("rotationalDegeneracy");

        if (!arguments.empty())
            Log<NumArgumentsError>::Error("rotationalDegeneracy");

        if (states.at(0)->type != rotational)
            Log<Message>::Error("Trying to assign rotational degeneracy to non-rotational state.");

        for (auto *state : states) {
            double J;

            if (!Parse::getValue(state->J, J))
                Log<Message>::Error("Non numerical rot level (" + state->J +
                                    ") when trying to assign rigid rotor energy.");

            state->statisticalWeight = 2 * J + 1;
        }
    }

    inline void constantValue(const std::vector<GasBase::StateBase*> &states,
                       double value, StatePropertyType type) {

//            if (arguments.size() != 1)
//                Log<NumArgumentsError>::Error("constantValue");
//
//            const double &value = arguments[0];

        for (auto *state : states) {
            setStateProperty(state, value, type);
        }
    }

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
}

#endif //LOKI_CPP_PROPERTYFUNCTIONS_H
