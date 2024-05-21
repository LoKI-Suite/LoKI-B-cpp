/** \file
 *
 *  Implementation of the Gas class.
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
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   May 2019
 */

#include "LoKI-B/Gas.h"
#include "LoKI-B/Log.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>

namespace loki
{

Gas::State::State(const StateEntry &entry, const Gas *gas, const State &parent)
    : m_gas(gas), m_parent(&parent), type(static_cast<StateType>(parent.type + 1)),
      charge(type == StateType::charge ? entry.m_charge : parent.charge),
      e(type == StateType::electronic ? entry.m_e : parent.e), v(type == StateType::vibrational ? entry.m_v : parent.v),
      J(type == StateType::rotational ? entry.m_J : parent.J), energy(-1), statisticalWeight(-1),
      m_population(0),
      m_delta(0)
{
    switch (type)
    {
    case StateType::electronic:
        if (e.empty())
            throw std::runtime_error("Electronic state not specified.");
        break;
    case StateType::vibrational:
        if (v.empty())
            throw std::runtime_error("Vibrational state not specified.");
        break;
    case StateType::rotational:
        if (J.empty())
            throw std::runtime_error("Rotational state not specified.");
        break;
    default:
        break;
    }
    /** \todo if multiple ionization levels are present, we set the population of the
     *        neutral state to 1 (the others remain at 0). EDIT: but only if the gas
     *        fraction is > 0. However, the gas fractions have not been loaded at this
     *        point, therefore this is now done in the "propagateFraction" function.
     */
    // if (m_gas->fraction != 0 && type == StateType::charge && charge.empty())
    // {
    //     setpopulation(1);
    // }
    assert(m_gas);
#if 0
    std::cout << "Created state: " << *this << ", type = " << type << std::endl;
    std::cout << " - entry: " << entry << std::endl;
    std::cout << " - parent: " << parent << ", type = " << parent.type << std::endl;
#endif
}

Gas::State::State(const Gas *gas)
    : m_gas(gas), m_parent(nullptr), type(StateType::root), charge(std::string{}), e(std::string{}),
      v(std::string{}), J(std::string{}), energy(-1), statisticalWeight(-1), m_population(0), m_delta(0)
{
    assert(m_gas);
#if 0
    std::cout << "State::State" << ", type = " << type << std::endl;
#endif
}

Gas::State::~State()
{
}

bool Gas::State::operator>=(const StateEntry &entry)
{
    if (type > entry.m_level)
        return false;

    if (charge == entry.m_charge)
    {
        if (type == StateType::charge)
            return true;
    }
    else
    {
        return false;
    }
    if (e == entry.m_e)
    {
        if (type == StateType::electronic)
            return true;
    }
    else
    {
        return false;
    }

    if (v == entry.m_v)
    {
        if (type == StateType::vibrational)
            return true;
    }
    else
    {
        return false;
    }

    if (J == entry.m_J)
    {
        if (type == StateType::rotational)
            return true;
    }
    else
    {
        return false;
    }

    return true;
}

std::ostream &operator<<(std::ostream &os, const Gas::State &state)
{
    os << state.gas().name();
#if 0
    // write N2+(...) instead of N2(+,...) ?
    /// \todo this does not work: when reading legacy input, the charge is smth. like '+', not a number.
    if (!state.charge.empty())
    {
        const int c=std::stoi(state.charge);
        os << std::string( std::abs(c), c > 0 ? '+' : '-');
    }
#endif

    // the electron is handled specially. The logic here is the same as
    // for in StateEntry's stream insertion operator.
    if (state.gas().name() == "e")
    {
        return os;
    }

    os << '(';

    if (!state.charge.empty())
    {
        os << state.charge << ',';
    }

    os << state.e;

    if (state.type >= StateType::vibrational)
        os << ",v=" << state.v;

    if (state.type == StateType::rotational)
        os << ",J=" << state.J;

    os << ')';

    return os;
}

void Gas::State::printChildren(std::ostream &os) const
{
    std::string space = "  ";

    for (uint8_t i = electronic; i < type; ++i)
        space.append("  ");

    for (auto state : children())
    {
        if (state != nullptr)
        {
            os << space << *state << std::endl;
            state->printChildren(os);
        }
    }
}

void Gas::State::checkPopulations() const
{
    if (children().empty())
        return;

    double totalPopulation = 0.;

    for (auto *state : children())
    {
        totalPopulation += state->population();
        state->checkPopulations();
    }

    if (population() == 0)
    {
        if (totalPopulation != 0)
            Log<ChildrenPopulationError>::Error(*this);
    }
    else
    {
        if (std::abs(totalPopulation - 1.) > 100. * std::numeric_limits<double>::epsilon())
            Log<ChildrenPopulationError>::Error(*this);
    }
}

void Gas::State::evaluateReducedDensities()
{
    /* \todo This old implementation appears to assign the (neutral) gas delta to the charged species.
     *       Setting that to zero for now...
     */
    if (type == StateType::root)
        m_delta = population() * gas().fraction;
    else
        m_delta = population() * parent()->delta();

    for (auto *state : children())
    {
        state->evaluateReducedDensities();
    }
}

Gas::State *Gas::State::find(const StateEntry &entry)
{
    auto it =
        std::find_if(children().begin(), children().end(), [&entry](State *child) { return *child >= entry; });
    return it == children().end() ? nullptr : *it;
}

Gas::Gas(const GasProperties& gasProps, std::string name)
    : m_root(new State(this)), m_name{name},
      mass{gasProps.get("mass",name)},
      harmonicFrequency{gasProps.get("harmonicFrequency",name,-1)},
      anharmonicFrequency{gasProps.get("anharmonicFrequency",name,-1)},
      rotationalConstant{gasProps.get("rotationalConstant",name,-1)},
      //electricDipoleMoment{gasProps.get("electricDipoleMoment",name,-1)},
      electricQuadrupoleMoment{gasProps.get("electricQuadrupoleMoment",name,-1)},
      //polarizability{gasProps.get("polarizability",name,-1)},
      fraction{gasProps.get("fraction",name,0)}
{
}

Gas::~Gas()
{
}

void Gas::print(std::ostream &os) const
{
    m_root->printChildren(os);
}

void Gas::propagateFraction()
{
    if (fraction > 0)
    {
        m_root->setPopulation(1.);
        auto it = std::find_if(m_root->children().begin(), m_root->children().end(),
                               [](auto *state) { return state->charge.empty(); });

        if (it != m_root->children().end())
        {
            Log<Message>::Notify("Set population of ", *(*it), " to one.");
            (*it)->setPopulation(1.);
        }
    }
}

void Gas::checkPopulations()
{
    m_root->checkPopulations();
}

void Gas::evaluateReducedDensities()
{
    m_root->evaluateReducedDensities();
}

Gas::State *Gas::findState(const StateEntry &entry)
{
    // Find charge state.

    auto state = m_root->find(entry);

    if (state == nullptr)
        return nullptr;

    if (entry.m_level == charge)
        return state;

    if (entry.m_e == "*" && entry.m_level == electronic)
    {
        if (state->children().empty())
            return nullptr;

        return state->children()[0];
    }

    // Find electronic state.

    state = state->find(entry);

    if (state == nullptr)
        return nullptr;

    if (entry.m_level == electronic)
        return state;

    if (entry.m_v == "*" && entry.m_level == vibrational)
    {
        if (state->children().empty())
            return nullptr;

        return state->children()[0];
    }

    // Find vibrational state.

    state = state->find(entry); // findState(state, entry);

    if (state == nullptr)
        return nullptr;

    if (entry.m_level == vibrational)
        return state;

    if (entry.m_J == "*" && entry.m_level == rotational)
    {
        if (state->children().empty())
            return nullptr;

        return state->children()[0];
    }

    // Find rotational state

    state = state->find(entry); // findState(state, entry);

    if (state == nullptr)
        return nullptr;

    if (entry.m_level == rotational)
        return state;

    return nullptr;
}

} // namespace loki
