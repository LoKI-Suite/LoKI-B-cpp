//
// Created by daan on 28-5-19.
//

#ifndef LOKI_CPP_PROPERTYFUNCTIONS_H
#define LOKI_CPP_PROPERTYFUNCTIONS_H

#include "Enumeration.h"
#include "Constant.h"
#include "Traits.h"
#include "Log.h"

#include <cmath>

// TODO: comment loki::PropertyFunctions

namespace loki::PropertyFunctions {
    using namespace Constant;

    template<typename TraitType>
    static void
    setStateProperty(typename Trait<TraitType>::State *state, const double &value, StatePropertyType type) {
        switch (type) {
            case StatePropertyType::energy:
                state->energy = value;
            case StatePropertyType::statisticalWeight:
                state->statisticalWeight = value;
            case StatePropertyType::population:
                state->population = value;
        }
    }

    template<typename TraitType>
    void boltzmannPopulation(std::vector<typename Trait<TraitType>::State *> &states,
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

    template<typename TraitType>
    void harmonicOscillatorEnergy(std::vector<typename Trait<TraitType>::State *> &states,
                                  const std::vector<double> &arguments, StatePropertyType type) {
        if (type != StatePropertyType::energy)
            Log<WrongPropertyError>::Error("harmonicOscillatorEnergy");

        if (!arguments.empty())
            Log<NumArgumentsError>::Error("harmonicOscillatorEnergy");

        if (states.at(0)->type != vibrational)
            Log<Message>::Error("Trying to assign harmonic oscillator energy to non-vibrational state.");

        if (states.at(0)->gas->harmonicFrequency < 0)
            Log<Message>::Error("Cannot find harmonicFrequency of the gas " + states.at(0)->gas->name +
                                " to evaluate state energies.");

        for (auto *state : states) {
            double vibLevel;

            if (!Parse::getValue(state->v, vibLevel))
                Log<Message>::Error("Non numerical vib level (" + state->v +
                                    ") when trying to assign harmonic oscillator energy.");

            state->energy = plankReducedInEv * state->gas->harmonicFrequency * (vibLevel + .5);
        }
    }

    template<typename TraitType>
    void rigidRotorEnergy(std::vector<typename Trait<TraitType>::State *> &states,
                          const std::vector<double> &arguments, StatePropertyType type) {
        if (type != StatePropertyType::energy)
            Log<WrongPropertyError>::Error("rigidRotorEnergy");

        if (!arguments.empty())
            Log<NumArgumentsError>::Error("rigidRotorEnergy");

        if (states.at(0)->type != rotational)
            Log<Message>::Error("Trying to assign rigid rotor energy to non-rotational state.");

        if (states.at(0)->gas->harmonicFrequency < 0)
            Log<Message>::Error("Cannot find rotationalConstant of the gas " + states.at(0)->gas->name +
                                " to evaluate state energies.");

        for (auto *state : states) {
            double rotLevel;

            if (!Parse::getValue(state->J, rotLevel))
                Log<Message>::Error("Non numerical rot level (" + state->J +
                                    ") when trying to assign rigid rotor energy.");

            state->energy = state->gas->rotationalConstant * rotLevel * (rotLevel + 1.);
        }
    }

    template<typename TraitType>
    void rotationalDegeneracy_N2(std::vector<typename Trait<TraitType>::State *> &states,
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

    template<typename TraitType>
    void constantValue(std::vector<typename Trait<TraitType>::State *> &states,
                       double value, StatePropertyType type) {

//            if (arguments.size() != 1)
//                Log<NumArgumentsError>::Error("constantValue");
//
//            const double &value = arguments[0];

        for (auto *state : states) {
            setStateProperty<TraitType>(state, value, type);
        }
    }

    template<typename TraitType>
    void callByName(const std::string &name, std::vector<typename Trait<TraitType>::State *> &states,
                    const std::vector<double> &arguments, StatePropertyType type) {

        if (name == "boltzmannPopulation")
            boltzmannPopulation<TraitType>(states, arguments, type);
        else if (name == "harmonicOscillatorEnergy")
            harmonicOscillatorEnergy<TraitType>(states, arguments, type);
        else if (name == "rigidRotorEnergy")
            rigidRotorEnergy<TraitType>(states, arguments, type);
        else if (name == "rotationalDegeneracy_N2")
            rotationalDegeneracy_N2<TraitType>(states, arguments, type);
        else
            Log<PropertyFunctionError>::Error(name);
    }
}

#endif //LOKI_CPP_PROPERTYFUNCTIONS_H
