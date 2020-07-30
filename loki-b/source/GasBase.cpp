#include "LoKI-B/GasBase.h"
#include "LoKI-B/Log.h"
#include <cassert>
#include <limits>
#include <algorithm>
#include <iostream>

namespace loki {

GasBase::StateBase::StateBase(const StateEntry &entry, StateType type, GasBase* gas, StateBase* parent)
: type(type), e(entry.e), v(entry.v), J(entry.J), charge(entry.charge), m_gas_base(gas), m_parent_base(parent),
  population(0), energy(-1), statisticalWeight(-1), density(0)
{
    assert(m_gas_base);
#if 0
    std::cout << "StateBase::StateBase from entry: " << std::endl
	<< " type = " << static_cast<int>(type) << std::endl
	<< " gas  = " << gas->name << std::endl
	<< " e    = " << e << std::endl
	<< " v    = " << v << std::endl
	<< " J    = " << J << std::endl
	<< " q    = " << charge << std::endl
        << std::endl;
#endif
}

GasBase::StateBase::StateBase(StateType type, std::string e, std::string v, std::string J, std::string charge, GasBase* gas)
: type(type), e(e), v(v), J(J), charge(charge), m_gas_base(gas), m_parent_base(nullptr),
  population(0), energy(-1), statisticalWeight(-1), density(0)
{
    assert(m_gas_base);
#if 0
    std::cout << "StateBase::StateBase from entry: " << std::endl
	<< " type = " << static_cast<int>(type) << std::endl
	<< " gas  = " << gas->name << std::endl
	<< " e    = " << e << std::endl
	<< " v    = " << v << std::endl
	<< " J    = " << J << std::endl
	<< " q    = " << charge << std::endl
        << std::endl;
#endif
}

GasBase::StateBase::~StateBase()
{
}

bool GasBase::StateBase::operator>=(const StateEntry &entry)
{
    {
        if (type > entry.level || charge != entry.charge)
            return false;

        if (e == entry.e)
        {
            if (type == Enumeration::electronic)
                return true;
        }
        else
        {
            return false;
        }

        if (v == entry.v)
        {
            if (type == Enumeration::vibrational)
                return true;
        }
        else
        {
            return false;
        }

        if (J == entry.J)
        {
            if (type == Enumeration::rotational)
                return true;
        }
        else
        {
            return false;
        }

        return true;
    }
}

std::ostream &operator<<(std::ostream &os, const GasBase::StateBase &state)
{
    os << state.gas_base().name << '(';

    if (!state.charge.empty())
        os << state.charge << ',';

    os << state.e;

    if (state.type >= Enumeration::vibrational)
        os << ",v=" << state.v;

    if (state.type == Enumeration::rotational)
        os << ",J=" << state.J;

    os << ')';

    return os;
}

void GasBase::StateBase::printChildren(std::ostream& os) const
{
    std::string space = "  ";

    for (uint8_t i = electronic; i < type; ++i)
        space.append("  ");

    for (auto state : m_children)
    {
        if (state != nullptr)
        {
            os << space << *state << std::endl;
            state->printChildren(os);
        }
    }
}

void GasBase::StateBase::checkPopulations() const
{
    if (m_children.empty())
        return;

    double totalPopulation = 0.;

    for (auto *state : m_children)
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

void GasBase::StateBase::evaluateDensity()
{
    if (type == electronic)
        density = population * gas_base().fraction;
    else
        density = population * parent_base()->density;

    for (auto *state : m_children)
        state->evaluateDensity();
}

GasBase::StateBase* GasBase::StateBase::find(const StateEntry &entry)
{
    auto it = std::find_if(m_children.begin(), m_children.end(),
                           [&entry](StateBase *child) { return *child >= entry; });
    return it == m_children.end() ? nullptr : *it;
}

GasBase:: GasBase(std::string name)
    : name{name},
    mass{-1},
    harmonicFrequency{-1},
    anharmonicFrequency{-1},
    rotationalConstant{-1},
    electricDipoleMoment{-1},
    electricQuadrupoleMoment{-1},
    polarizability{-1},
    fraction{0}
{
}

void GasBase::print(std::ostream& os) const
{
    for (const auto& s : stateBaseTree)
    {
        os << *s << std::endl;
        s->printChildren(os);
    }
    for (const auto& s : ionicBaseStates)
    {
        os << *s << std::endl;
        s->printChildren(os);
    }
}

void GasBase::checkPopulations()
{
    double totalPopulation = 0.;

    for (auto *state : stateBaseTree) {
        totalPopulation += state->population;
        state->checkPopulations();
    }
    for (auto *state : ionicBaseStates) {
        totalPopulation += state->population;
        state->checkPopulations();
    }

    if (fraction == 0) {
        if (totalPopulation != 0)
            Log<ZeroFractionPopulationError>::Error(name);
    } else if (std::abs(totalPopulation - 1.) > 10. * std::numeric_limits<double>::epsilon()) {
        Log<ChildrenPopulationError>::Error(name);
    }
}

void GasBase::evaluateStateDensities()
{
    for (auto *state : stateBaseTree)
        state->evaluateDensity();
    for (auto *state : ionicBaseStates)
        state->evaluateDensity();
}

GasBase::StateBase* GasBase::find(const StateEntry& entry)
{
    auto &childStates = (entry.charge.empty() ? stateBaseTree : ionicBaseStates);

    auto it = std::find_if(childStates.begin(), childStates.end(),
                           [&entry](StateBase* state) {
                               return *state >= entry;
                           });

    if (it == childStates.end()) {
        return nullptr;
    }

    return *it;
}

GasBase::StateBase* GasBase::findState(const StateEntry &entry)
{
    if (entry.e == "*" && entry.level == electronic) {
        auto &states = entry.charge.empty() ? stateBaseTree : ionicBaseStates;

        if (states.empty())
            return nullptr;

        return states[0];
    }

    // Find electronic state.

    auto * state = find(entry);

    if (state == nullptr)
        return nullptr;

    if (entry.level == electronic)
        return state;

    if (entry.v == "*" && entry.level == vibrational) {
        if (state->m_children.empty())
            return nullptr;

        return state->m_children[0];
    }

    // Find vibrational state.

    state = state->find(entry);// findState(state, entry);

    if (state == nullptr)
        return nullptr;

    if (entry.level == vibrational)
        return state;

    if (entry.J == "*" && entry.level == rotational) {
        if (state->m_children.empty())
            return nullptr;

        return state->m_children[0];
    }

    // Find rotational state

    state = state->find(entry);//findState(state, entry);

    if (state == nullptr)
        return nullptr;

    if (entry.level == rotational)
        return state;

    return nullptr;
}

} // namespace loki
