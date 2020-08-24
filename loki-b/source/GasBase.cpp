#include "LoKI-B/GasBase.h"
#include "LoKI-B/Log.h"
#include <cassert>
#include <limits>
#include <algorithm>
#include <iostream>

namespace loki {

GasBase::StateBase::StateBase(const StateEntry &entry, GasBase* gas, StateBase& parent)
    : m_gas_base(gas),
    m_parent_base(&parent),
    type(static_cast<StateType>(parent.type + 1)),
    charge(type==StateType::charge ? entry.charge : parent.charge),
    e(type==StateType::electronic ? entry.e : parent.e),
    v(type==StateType::vibrational ? entry.v : parent.v),
    J(type==StateType::rotational ? entry.J : parent.J),
    population(0),
    energy(-1),
    statisticalWeight(-1),
    density(0)
{
    switch (type)
    {
        case StateType::electronic:
            if (e.empty()) throw std::runtime_error("Electronic state not specified.");
        break;
        case StateType::vibrational:
            if (v.empty()) throw std::runtime_error("Vibrational state not specified.");
        break;
        case StateType::rotational:
            if (J.empty()) throw std::runtime_error("Rotational state not specified.");
        break;
    }
    /** \todo if multiple ionization levels are present, we set the population of the
     *        neutral state to 1 (the others remain at 0).
     */
    if (type==StateType::charge && charge.empty())
    {
        population=1;
    }
    assert(m_gas_base);
#if 0
    std::cout << "Created state: " << *this << ", type = " << type << std::endl;
    std::cout << " - entry: " << entry << std::endl;
    std::cout << " - parent: " << parent << ", type = " << parent.type << std::endl;
#endif
}

GasBase::StateBase::StateBase(GasBase* gas)
    : m_gas_base(gas),
    m_parent_base(nullptr),
    type(StateType::root),
    charge(std::string{}),
    e(std::string{}),
    v(std::string{}),
    J(std::string{}),
    population(1),
    energy(-1),
    statisticalWeight(-1),
    density(0)
{
    assert(m_gas_base);
#if 0
    std::cout << "StateBase::StateBase" << ", type = " << type << std::endl;
#endif
}

GasBase::StateBase::~StateBase()
{
}

bool GasBase::StateBase::operator>=(const StateEntry &entry)
{
    if (type > entry.level)
        return false;

    if (charge == entry.charge)
    {
        if (type == StateType::charge)
            return true;
    }
    else
    {
        return false;
    }
    if (e == entry.e)
    {
        if (type == StateType::electronic)
            return true;
    }
    else
    {
        return false;
    }

    if (v == entry.v)
    {
        if (type == StateType::vibrational)
            return true;
    }
    else
    {
        return false;
    }

    if (J == entry.J)
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

std::ostream &operator<<(std::ostream &os, const GasBase::StateBase &state)
{
    os << state.gas().name;
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
    if (state.gas().name=="e")
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
        if (std::abs(totalPopulation - 1.) > 100. * std::numeric_limits<double>::epsilon())
            Log<ChildrenPopulationError>::Error(*this);
    }
}

void GasBase::StateBase::evaluateDensity()
{
    /* \todo This old implementation appears to assign the (neutral) gas density to the charged species.
     *       Setting that to zero for now...
     */
    if (type == StateType::root)
        density = population * gas().fraction;
    else
        density = population * parent()->density;

    for (auto *state : m_children)
        state->evaluateDensity();
}

GasBase::StateBase* GasBase::StateBase::find(const StateEntry &entry)
{
    auto it = std::find_if(m_children.begin(), m_children.end(),
                           [&entry](StateBase *child) { return *child >= entry; });
    return it == m_children.end() ? nullptr : *it;
}

GasBase::GasBase(std::string name)
    : m_root(new State(this)),
    name{name},
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

GasBase::~GasBase()
{
}

void GasBase::print(std::ostream& os) const
{
    m_root->printChildren(os);
}

void GasBase::checkPopulations()
{
    m_root->checkPopulations();
}

void GasBase::evaluateStateDensities()
{
    m_root->evaluateDensity();
}

GasBase::StateBase* GasBase::findState(const StateEntry &entry)
{
    // Find charge state.

    auto state = m_root->find(entry);

    if (state == nullptr)
        return nullptr;

    if (entry.level == charge)
        return state;

    if (entry.e == "*" && entry.level == electronic) {
        if (state->m_children.empty())
            return nullptr;

        return state->m_children[0];
    }

    // Find electronic state.

    state = state->find(entry);

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
