/** \file
 *
 *  Declaration of functions that describe gas and state properties.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \author Daan Boer
 *  \date   28 May 2019
 */

#ifndef LOKI_CPP_PROPERTYFUNCTIONS_H
#define LOKI_CPP_PROPERTYFUNCTIONS_H

#include "LoKI-B/Constant.h"
#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Gas.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/Enumeration.h"

#include <cmath>

/** This namespace hosts all property functions. Each property function has
 *  the same structure.
 *
 *  The arguments are as follows.
 *
 *  1. A vector of pointers to the states (State) of which to set the
       property.
 *  2. A vector of doubles representing the arguments of the property function.
 *     As an example, in the case of boltzmannPopulation this vector contains a
 *     single argument, the temperature.
 *  3. A StatePropertyType indicating whether the function is supposed to set
 *     either the energy, population or statistical weight of the states.
 *
 *  The user is responsible to provide any checks on the number of arguments
 *  provided, the validity of the passed StateProperty etc. Furthermore,
 *  when adding a function, the user should also add it to the callByName
 *  function, which links the function to a string containing its name. When
 *  this is done, it is possible to address the function in a setup file to
 *  set state properties.
 */
namespace loki::PropertyFunctions
{

/** Accepts a pointer to a state, a double value and a StatePropertyType. Based on the value of the
 *  StatePropertyType, it will assign the value to either the population, statisticalWeight or
 *  population member variables of the supplied state object.
 */
inline void setStateProperty(Gas::State *state, const double &value, StatePropertyType type)
{
    switch (type)
    {
    case StatePropertyType::energy:
        state->energy = value;
        break;
    case StatePropertyType::statisticalWeight:
        state->statisticalWeight = value;
        break;
    case StatePropertyType::population:
        state->setPopulation(value);
        break;
    }
}

inline void boltzmannPopulation(const std::vector<Gas::State *> &states, const std::vector<double> &arguments,
                                StatePropertyType type)
{
    if (type != StatePropertyType::population)
        Log<WrongPropertyError>::Error("boltzmannPopulation");

    if (arguments.size() != 1)
        Log<NumArgumentsError>::Error("boltzmannPopulation");

    const double &temp = arguments[0];
    double norm = 0.;

    for (auto *state : states)
    {
        state->setPopulation(state->statisticalWeight * std::exp(-state->energy/(Constant::kBeV * temp)));
        norm += state->population();
    }
    for (auto *state : states)
    {
        state->setPopulation(state->population()/norm);
    }
}

inline void boltzmannPopulationVibrationalCutoff(const std::vector<Gas::State *> &states, const std::vector<double> &arguments,
                                                 StatePropertyType type)
{
    if (type != StatePropertyType::population)
        Log<WrongPropertyError>::Error("boltzmannPopulationVibrationalCutoff");

    if (arguments.size() != 2)
        Log<NumArgumentsError>::Error("boltzmannPopulationVibrationalCutoff");

    const double &temp = arguments[0];
    const double &vMax = arguments[1];
    double norm = 0.;

    for (auto *state : states)
    {
        double vibLevel;

        if (state->type != StateType::vibrational)
            Log<Message>::Error("Trying to assign a Boltzmann population with vibrational cutoff to a non-vibrational state.");

        if (!Parse::getValue(state->v, vibLevel))
            Log<Message>::Error("Non numerical vib level (" + state->v +
                                ") when trying to assign Boltzmann population with vibrational cutoff.");

        if (vibLevel > vMax) {
            state->setPopulation(0);
        } else {
            state->setPopulation(state->statisticalWeight * std::exp(-state->energy/(Constant::kBeV * temp)));
            norm += state->population();
        }
    }
    for (auto *state : states)
    {
        state->setPopulation(state->population()/norm);
    }
}

inline void boltzmannPopulationRotationalCutoff(const std::vector<Gas::State *> &states, const std::vector<double> &arguments,
                                                 StatePropertyType type)
{
    if (type != StatePropertyType::population)
        Log<WrongPropertyError>::Error("boltzmannPopulationRotationalCutoff");

    if (arguments.size() != 2)
        Log<NumArgumentsError>::Error("boltzmannPopulationRotationalCutoff");

    const double &temp = arguments[0];
    const double &jMax = arguments[1];
    double norm = 0.;

    for (auto *state : states)
    {
        double rotLevel;

        if (state->type != StateType::rotational)
            Log<Message>::Error("Trying to assign a Boltzmann population with rotational cutoff to a non-rotational state.");

        if (!Parse::getValue(state->J, rotLevel))
            Log<Message>::Error("Non numerical rot level (" + state->v +
                                ") when trying to assign Boltzmann population with rotational cutoff.");

        if (rotLevel > jMax) {
            state->setPopulation(0);
        } else {
            state->setPopulation(state->statisticalWeight * std::exp(-state->energy/(Constant::kBeV * temp)));
            norm += state->population();
        }
    }
    for (auto *state : states)
    {
        state->setPopulation(state->population()/norm);
    }
}

inline void treanorPopulation(const std::vector<Gas::State *> &states, const std::vector<double> &arguments,
                              StatePropertyType type)
{
    if (type != StatePropertyType::population)
        Log<WrongPropertyError>::Error("treanorPopulation");

    if (arguments.size() != 2)
        Log<NumArgumentsError>::Error("treanorPopulation");

    const double &tempZero = arguments[0];
    const double &tempOne = arguments[1];
    double eZero = -1;
    double eOne = -1;
    double norm = 0.;

    for (const auto *state : states) {
        if (state->type != StateType::vibrational)
            Log<Message>::Error("Trying to assign a Treanor population to non-vibrational state.");

        if (state->v == "0") {
            eZero = state->energy;
        } else if (state->v == "1") {
            eOne = state->energy;
        }
    }

    if (eZero == -1 || eOne == -1)
        Log<Message>::Error("Unable to find E_0 or E_1 to apply Treanor population.");

    for (auto *state : states)
    {
        double vibLevel;

        if (!Parse::getValue(state->v, vibLevel))
            Log<Message>::Error("Non numerical vib level (" + state->v +
                                ") when trying to assign Treanor population.");

        const double energy = state->energy;
        const double g = state->statisticalWeight;

        if (energy == -1)
            Log<Message>::Error("Unable to find energy of state while applying Treanor population.");
        if (g == -1)
            Log<Message>::Error("Unable to find statistical weight of state while applying Treanor population.");

        state->setPopulation(g * std::exp(-(vibLevel * (eOne - eZero) * (1/tempOne - 1/tempZero) + (energy - eZero) / tempZero) / Constant::kBeV));
        norm += state->population();
    }
    for (auto *state : states)
    {
        state->setPopulation(state->population()/norm);
    }
}

inline void treanorGordietsPopulation(const std::vector<Gas::State *> &states, const std::vector<double> &arguments,
                                      StatePropertyType type)
{
    if (type != StatePropertyType::population)
        Log<WrongPropertyError>::Error("treanorGordietsPopulation");

    if (arguments.size() != 2)
        Log<NumArgumentsError>::Error("treanorGordietsPopulation");

    const double &tempZero = arguments[0];
    const double &tempOne = arguments[1];
    double eZero = -1;
    double eOne = -1;
    double norm = 0.;

    for (const auto *state : states) {
        if (state->type != StateType::vibrational)
            Log<Message>::Error("Trying to assign a Treanor-Gordiets population to non-vibrational state.");

        if (state->v == "0") {
            eZero = state->energy;
        } else if (state->v == "1") {
            eOne = state->energy;
        }
    }

    if (eZero == -1 || eOne == -1)
        Log<Message>::Error("Unable to find E_0 or E_1 to apply Treanor-Gordiets population.");

    if (states.at(0)->gas().anharmonicFrequency < 0)
        Log<Message>::Error("Cannot find anharmonicFrequency of the gas " + states.at(0)->gas().name() +
                            " to evaluate state energies.");

    const double vLimit = std::floor(0.5 * (1 + (eOne - eZero) * tempZero / (Constant::plankReducedInEv * states.at(0)->gas().anharmonicFrequency * tempOne)));
    double vLimitPop = 0;

    const auto computePopulation = [eZero, eOne, tempZero, tempOne](double v, double g, double energy){
        return g * std::exp(-(v * (eOne - eZero) * (1/tempOne - 1/tempZero) + (energy - eZero) / tempZero) / Constant::kBeV);
    };

    for (auto *state : states)
    {
        double vibLevel;

        if (!Parse::getValue(state->v, vibLevel))
            Log<Message>::Error("Non numerical vib level (" + state->v +
                                ") when trying to assign Treanor-Gordiets population.");

        const double energy = state->energy;
        const double g = state->statisticalWeight;

        if (energy == -1)
            Log<Message>::Error("Unable to find energy of state while applying Treanor-Gordiets population.");
        if (g == -1)
            Log<Message>::Error("Unable to find statistical weight of state while applying Treanor-Gordiets population.");

        if (vibLevel < vLimit) {
          state->setPopulation(computePopulation(vibLevel, g, energy));
        } else if (vibLevel == vLimit) {
          vLimitPop = computePopulation(vibLevel, g, energy);
          state->setPopulation(vLimitPop);
        } else {
          state->setPopulation(vLimitPop * vLimit / vibLevel);
        }
        norm += state->population();
    }
    for (auto *state : states)
    {
        state->setPopulation(state->population()/norm);
    }
}

inline void harmonicOscillatorEnergy(const std::vector<Gas::State *> &states,
                                     const std::vector<double> &arguments, StatePropertyType type)
{
    if (type != StatePropertyType::energy)
        Log<WrongPropertyError>::Error("harmonicOscillatorEnergy");

    if (!arguments.empty())
        Log<NumArgumentsError>::Error("harmonicOscillatorEnergy");

    if (states.at(0)->type != vibrational)
        Log<Message>::Error("Trying to assign harmonic oscillator energy to non-vibrational state.");

    if (states.at(0)->gas().harmonicFrequency < 0)
        Log<Message>::Error("Cannot find harmonicFrequency of the gas " + states.at(0)->gas().name() +
                            " to evaluate state energies.");

    for (auto *state : states)
    {
        double vibLevel;

        if (!Parse::getValue(state->v, vibLevel))
            Log<Message>::Error("Non numerical vib level (" + state->v +
                                ") when trying to assign harmonic oscillator energy.");

        state->energy = Constant::plankReducedInEv * state->gas().harmonicFrequency * (vibLevel + .5);
    }
}

inline void morseOscillatorEnergy(const std::vector<Gas::State *> &states,
                                  const std::vector<double> &arguments, StatePropertyType type)
{
    if (type != StatePropertyType::energy)
        Log<WrongPropertyError>::Error("morseOscillatorEnergy");

    if (!arguments.empty())
        Log<NumArgumentsError>::Error("morseOscillatorEnergy");

    if (states.at(0)->type != vibrational)
        Log<Message>::Error("Trying to assign Morse oscillator energy to non-vibrational state.");

    if (states.at(0)->gas().harmonicFrequency < 0)
        Log<Message>::Error("Cannot find harmonicFrequency of the gas " + states.at(0)->gas().name() +
                            " to evaluate state energies.");

    if (states.at(0)->gas().anharmonicFrequency < 0)
        Log<Message>::Error("Cannot find anharmonicFrequency of the gas " + states.at(0)->gas().name() +
                            " to evaluate state energies.");

    for (auto *state : states)
    {
        double vibLevel;

        if (!Parse::getValue(state->v, vibLevel))
            Log<Message>::Error("Non numerical vib level (" + state->v +
                                ") when trying to assign Morse oscillator energy.");

        state->energy = Constant::plankReducedInEv * (state->gas().harmonicFrequency * (vibLevel + .5) - state->gas().anharmonicFrequency * std::pow(vibLevel + 0.5, 2));
    }
}

inline void rigidRotorEnergy(const std::vector<Gas::State *> &states, const std::vector<double> &arguments,
                             StatePropertyType type)
{
    if (type != StatePropertyType::energy)
        Log<WrongPropertyError>::Error("rigidRotorEnergy");

    if (!arguments.empty())
        Log<NumArgumentsError>::Error("rigidRotorEnergy");

    if (states.at(0)->type != rotational)
        Log<Message>::Error("Trying to assign rigid rotor energy to non-rotational state.");

    if (states.at(0)->gas().rotationalConstant < 0)
        Log<Message>::Error("Cannot find rotationalConstant of the gas " + states.at(0)->gas().name() +
                            " to evaluate state energies.");

    for (auto *state : states)
    {
        double rotLevel;

        if (!Parse::getValue(state->J, rotLevel))
            Log<Message>::Error("Non numerical rot level (" + state->J + ") when trying to assign rigid rotor energy.");

        state->energy = state->gas().rotationalConstant * rotLevel * (rotLevel + 1.);
    }
}

inline void rotationalDegeneracy_H2(const std::vector<Gas::State *> &states,
                                    const std::vector<double> &arguments, StatePropertyType type)
{
    if (type != StatePropertyType::statisticalWeight)
        Log<WrongPropertyError>::Error("rotationalDegeneracy_H2");

    if (!arguments.empty())
        Log<NumArgumentsError>::Error("rotationalDegeneracy_H2");

    if (states.at(0)->type != rotational)
        Log<Message>::Error("Trying to assign rotational degeneracy to non-rotational state.");

    for (auto *state : states)
    {
        double J;

        if (!Parse::getValue(state->J, J))
            Log<Message>::Error("Non numerical rot level (" + state->J + ") when trying to assign rigid rotor energy.");

        // See eq. 72 in \cite Manual_2_2_0. This expression is originally taken from 
        // \cite Ridenti2015
        state->statisticalWeight = (2 - std::pow(-1, J)) * (2 * J + 1);
    }
}

inline void rotationalDegeneracy_N2(const std::vector<Gas::State *> &states,
                                    const std::vector<double> &arguments, StatePropertyType type)
{
    if (type != StatePropertyType::statisticalWeight)
        Log<WrongPropertyError>::Error("rotationalDegeneracy_N2");

    if (!arguments.empty())
        Log<NumArgumentsError>::Error("rotationalDegeneracy_N2");

    if (states.at(0)->type != rotational)
        Log<Message>::Error("Trying to assign rotational degeneracy to non-rotational state.");

    for (auto *state : states)
    {
        double J;

        if (!Parse::getValue(state->J, J))
            Log<Message>::Error("Non numerical rot level (" + state->J + ") when trying to assign rigid rotor energy.");

        // See eq. 71 in \cite Manual_2_2_0. This expression is originally taken from 
        // \cite Ridenti2015
        state->statisticalWeight = 3 * (1. + .5 * (1. + std::pow(-1., J))) * (2. * J + 1.);
    }
}

inline void rotationalDegeneracy(const std::vector<Gas::State *> &states, const std::vector<double> &arguments,
                                 StatePropertyType type)
{
    if (type != StatePropertyType::statisticalWeight)
        Log<WrongPropertyError>::Error("rotationalDegeneracy");

    if (!arguments.empty())
        Log<NumArgumentsError>::Error("rotationalDegeneracy");

    if (states.at(0)->type != rotational)
        Log<Message>::Error("Trying to assign rotational degeneracy to non-rotational state.");

    for (auto *state : states)
    {
        double J;

        if (!Parse::getValue(state->J, J))
            Log<Message>::Error("Non numerical rot level (" + state->J + ") when trying to assign rigid rotor energy.");

        // See eq. 70 in \cite Manual_2_2_0. This expression is originally taken from 
        // \cite Ridenti2015
        state->statisticalWeight = 2 * J + 1;
    }
}

inline void constantValue(const std::vector<Gas::State *> &states, double value, StatePropertyType type)
{

    //            if (arguments.size() != 1)
    //                Log<NumArgumentsError>::Error("constantValue");
    //
    //            const double &value = arguments[0];

    for (auto *state : states)
    {
        setStateProperty(state, value, type);
    }
}

inline void callByName(const std::string &name, const std::vector<Gas::State *> &states,
                       const std::vector<double> &arguments, StatePropertyType type)
{
    if (name == "boltzmannPopulation")
        boltzmannPopulation(states, arguments, type);
    else if (name == "boltzmannPopulationVibrationalCutoff")
        boltzmannPopulationVibrationalCutoff(states, arguments, type);
    else if (name == "boltzmannPopulationRotationalCutoff")
        boltzmannPopulationRotationalCutoff(states, arguments, type);
    else if (name == "treanorPopulation")
        treanorPopulation(states, arguments, type);
    else if (name == "treanorGordietsPopulation")
        treanorGordietsPopulation(states, arguments, type);
    else if (name == "harmonicOscillatorEnergy")
        harmonicOscillatorEnergy(states, arguments, type);
    else if (name == "morseOscillatorEnergy")
        morseOscillatorEnergy(states, arguments, type);
    else if (name == "rigidRotorEnergy")
        rigidRotorEnergy(states, arguments, type);
    else if (name == "rotationalDegeneracy_H2")
        rotationalDegeneracy_H2(states, arguments, type);
    else if (name == "rotationalDegeneracy_N2")
        rotationalDegeneracy_N2(states, arguments, type);
    else if (name == "rotationalDegeneracy")
        rotationalDegeneracy(states, arguments, type);
    else
        Log<PropertyFunctionError>::Error(name);
}
} // namespace loki::PropertyFunctions

#endif // LOKI_CPP_PROPERTYFUNCTIONS_H
