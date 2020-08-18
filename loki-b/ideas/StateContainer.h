#ifndef LOKI_CPP_STATECONTAINER_H
#define LOKI_CPP_STATECONTAINER_H

#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include <limits>
#include <type_traits>

#include "LoKI-B/json.h"

namespace loki {
namespace experimental {

template <class ChildT>
using ChildContainer = std::vector<std::unique_ptr<ChildT>>;

/* Template arguments:
 *
 *  - Node data (other than parent, children etc.)
 *  - The Parent, Current and Child types. StateT must be derived from
 *    State<StateInfoT,ParentT,StateT,ChildT>. This is checked with a
 *    static_assert in the constructor.
 *
 *  NOTE the specialization for ChildT=void, which terminates the state sequence.
 */
template <typename StateInfoT, typename ParentT, typename StateT, typename ChildT>
class State : public StateInfoT
{
public:
    using StateInfo = StateInfoT;
    using Parent = ParentT;
    using Child = ChildT;
    using ChildContainer = loki::experimental::ChildContainer<Child>;

    using size_type = typename ChildContainer::size_type;
    using iterator = typename ChildContainer::iterator;
    using const_iterator = typename ChildContainer::const_iterator;

    State(const Parent& parent, const std::string& name)
    : StateInfo(name), m_parent(parent)
    {
        static_assert(std::is_base_of_v<State,StateT>);
    }

    const Parent& parent() const { return m_parent; }

    const ChildContainer& children() const { return m_children; }
    /** \todo We could just expose children(), but probably it
     *  is better to use indicet_iterator to present the children
     *  by reference instead of pointer. That will also be more
     *  const-safe.
     */
    bool empty() const { return m_children.empty(); }
    size_type size() const { return m_children.size(); }
    iterator begin() { return m_children.begin(); }
    iterator end() { return m_children.end(); }
    const_iterator begin() const { return m_children.begin(); }
    const_iterator end() const { return m_children.end(); }

    Child* ensure_state(const std::string& info)
    {
        Child* c = find(info);
        return c ? c : add_state(info);
    }
    template <class ...Args>
    auto ensure_state(const std::string& info, Args... args)
    {
        // create child state 'info' (when necessary),
        // call ensure_state on that object with the remaining args.
        return ensure_state(info)->ensure_state(args...);
    }
    Child* find(const std::string& info) const
    {
        for (const auto& c : m_children)
        {
            /// \todo This will call StateInfo::operator==(info), which is a bit too subtle.
            if (*c==info)
            {
                return c.get();
            }
        }
        return nullptr;
    }
private:
    Child* add_state(const std::string& info)
    {
        return m_children.emplace_back(new Child(*static_cast<StateT*>(this),info)).get();
    }
    const Parent& m_parent;
    ChildContainer m_children;
};

/* Specialization for ChildT=void. This does not have/support child states.
 */
template <typename StateInfoT, typename ParentT, typename StateT>
class State<StateInfoT,ParentT,StateT,void> : public StateInfoT
{
public:
    using StateInfo = StateInfoT;
    using Parent = ParentT;
    using Child = void;

    State(const Parent& parent, const std::string& name)
    : StateInfo(name), m_parent(parent)
    {
        static_assert(std::is_base_of_v<State,StateT>);
    }
    const Parent& parent() const { return m_parent; }
private:
    const Parent& m_parent;
};

/** \todo See if we can be smarter about the various overloads, there
 *  are a lot if also Operator can be a const/non-const reference.
 *  Pass all arguments by value and let the caller use std::ref, std::cref?
 */

template <class StateT, class Operation, std::enable_if_t<!std::is_same_v<typename StateT::Child,void>>* =nullptr>
void apply(const StateT& state, const Operation& operation, bool recurse)
{
    operation(state);
    if (recurse)
    {
        for (const auto& c : state)
        {
            apply(*c,operation,recurse);
        }
    }
}

// version for Child==void
template <class StateT, class Operation, std::enable_if_t<std::is_same_v<typename StateT::Child,void>>* =nullptr>
void apply(const StateT& state, const Operation& operation, bool recurse)
{
    operation(state);
}

template <class StateT, class Operation, std::enable_if_t<!std::is_same_v<typename StateT::Child,void>>* =nullptr>
void apply(StateT& state, const Operation& operation, bool recurse)
{
    if (recurse)
    {
        for (const auto& c : state)
        {
            apply(*c,operation,recurse);
        }
    }
}

// version for Child==void
template <class StateT, class Operation, std::enable_if_t<std::is_same_v<typename StateT::Child,void>>* =nullptr>
void apply(StateT& state, const Operation& operation, bool recurse)
{
    operation(state);
}

// common State data. Use accessors so we can implement some ops later using member function pointers.

class StateData
{
public:
    StateData()
    : m_population{0},
    m_energy{-1},
    m_statisticalWeight{-1},
    m_density{0}
    {}

    double population() const { return m_population; }
    double energy() const { return m_energy; }
    double statisticalWeight() const { return m_statisticalWeight; }
    double density() const { return m_density; }

    double& population() { return m_population; }
    double& energy() { return m_energy; }
    double& statisticalWeight() { return m_statisticalWeight; }
    double& density() { return m_density; }

private:
    double m_population;
    double m_energy;
    double m_statisticalWeight;
    double m_density;
};

// info for the various states: Root, Charge, Electronic, Vibrational, Rotational.

class RootInfo : public StateData
{
public:
    RootInfo(const std::string& gasname)
     : m_gasname(gasname)
    {
        /// \todo think harder about population initialization.
        population()=1.0;
    }
    bool operator==(const std::string& info) const
    {
        return m_gasname==info;
    }
    std::string str() const { return m_gasname; }
    const std::string& gasName() const { return m_gasname; }
private:
    const std::string m_gasname;
};

class ChargeInfo : public StateData
{
public:
    ChargeInfo(const json_type& info) : m_q(info.at("charge").get<std::string>()) {};
    ChargeInfo(const std::string& info) : m_q(info) {};
    bool operator==(const std::string& info) const
    {
        return m_q==info;
    }
    std::string str() const { return "charge="+m_q; }
    const std::string& q() const { return m_q; }
private:
    const std::string m_q;
};

class ElectronicInfo : public StateData
{
public:
    ElectronicInfo(const json_type& info) : m_e(info.at("e").get<std::string>()) {};
    ElectronicInfo(const std::string& info) : m_e(info) {};
    bool operator==(const std::string& info) const
    {
        return m_e==info;
    }
    std::string str() const { return "e="+m_e; }
    const std::string& e() const { return m_e; }
private:
    const std::string m_e;
};

class VibrationalInfo : public StateData
{
public:
    VibrationalInfo(const json_type& info) : m_v(info.at("v").get<std::string>()) {};
    VibrationalInfo(const std::string& info) : m_v(info) {};
    bool operator==(const std::string& info) const
    {
        return m_v==info;
    }
    std::string str() const { return "v="+m_v; }
    const std::string& v() const { return m_v; }
private:
    std::string m_v;
};

class RotationalInfo : public StateData
{
public:
    RotationalInfo(const json_type& info) : m_J(info.at("J").get<std::string>()) {};
    RotationalInfo(const std::string& info) : m_J(info) {};
    bool operator==(const std::string& info) const
    {
        return m_J==info;
    }
    std::string str() const { return "J="+m_J; }
    const std::string& J() const { return m_J; }
private:
    const std::string m_J;
};



class Gas;
class StateRoot;
class StateCharge;
class StateElectronic;
class StateVibrational;
class StateRotational;

// here we assemble the various pieces of the chain.
// Root has the Gas as parent (todo: reconsider that), StateCharge as parent.

class StateRoot : public State<RootInfo,Gas,StateRoot,StateCharge>
{
public:
    StateRoot(Gas& gas);
    const Gas& getGas() const { return parent(); }
    /** \todo see if we can provide level() in a better way without having this member in every State type below
     *  Same for getGas(). Make this a free function template?
     */
    static constexpr unsigned level() { return 0; }
};

class StateCharge : public State<ChargeInfo,StateRoot,StateCharge,StateElectronic>
{
public:
    using State::State;
    const Gas& getGas() const { return parent().getGas(); }
    static constexpr unsigned level() { return Parent::level()+1; }
};

class StateElectronic : public State<ElectronicInfo,StateCharge,StateElectronic,StateVibrational>
{
public:
    using State::State;
    const Gas& getGas() const { return parent().getGas(); }
    static constexpr unsigned level() { return Parent::level()+1; }
};

class StateVibrational : public State<VibrationalInfo,StateElectronic,StateVibrational,StateRotational>
{
public:
    using State::State;
    const Gas& getGas() const { return parent().getGas(); }
    static constexpr unsigned level() { return Parent::level()+1; }
};

class StateRotational : public State<RotationalInfo,StateVibrational,StateRotational,void>
{
public:
    using State::State;
    const Gas& getGas() const { return parent().getGas(); }
    static constexpr unsigned level() { return Parent::level()+1; }
};

class Gas
{
public:
    explicit Gas(const std::string& name);
    /** arguments: C [, E [, V [, J ] ] ],
     * where C, E, V and J are info-strings for the charge, electronic,
     * vibrational and rotational aspects of the state. As an example,
     * ensure_state("0","X") creates electronic state N2(X), while the four-argument
     * invocation ensure_state("0","X","0","2") creates the neutral state
     * X(v=2,J=2).
     */
    template <class ...Args>
    auto ensure_state(Args... args)
    {
        return m_states->ensure_state(args...);
    }
    void printStates(std::ostream& os) const;
    const std::string& name() const { return m_name; }

    void checkPopulations() const;
    void evaluateStateDensities() const;

    double harmonicFrequency;
    double rotationalConstant;
    double fraction;

    const StateRoot& states() const { return *m_states; }
    StateRoot& states() { return *m_states; }
private:
    const std::string m_name;
    std::unique_ptr<StateRoot> m_states;
};


// implementation:

inline StateRoot::StateRoot(Gas& gas)
    : State<RootInfo,Gas,StateRoot,StateCharge>(gas,gas.name())
{
}

inline Gas::Gas(const std::string& name)
    : m_name(name),
    m_states(new StateRoot(*this)), harmonicFrequency(-1), rotationalConstant(-1), fraction(1)
{
}

inline void Gas::printStates(std::ostream& os) const
{
    apply(states(), [&](const auto& st) { os << std::string(st.level(),' ') << st.str() << std::endl; }, true );
}

struct PopulationChecker
{
    void operator()(const StateRotational& state) const
    {
        // nothing, does not have children
    }
    template <class StateT>
    void operator()(const StateT& state) const
    {
        if (state.empty())
        {
            return;
        }
        double totalPopulation = 0.;
        for (auto& c : state)
        {
            totalPopulation += c->population();
        }
        if (std::abs(totalPopulation - 1.) > 10. * std::numeric_limits<double>::epsilon())
        {
            throw std::runtime_error("Populations of children of '" + state.str() + "' do not add up to unity.");
        }
    }
};

inline void Gas::checkPopulations() const
{
    std::cout << "Checking populations of gas '" << m_name << std::endl;
    apply(states(),PopulationChecker{}, true);
}

struct DensityEvaluator
{
    void operator()(StateRoot& state) const
    {
        state.density() = state.population() * state.getGas().fraction;
    }
    template <class StateT>
    void operator()(StateT& state) const
    {
        state.density() = state.population() * state.parent().density();
    }
};

inline void Gas::evaluateStateDensities() const
{
    std::cout << "Evaluating densities of gas '" << m_name << std::endl;
    apply(*m_states,DensityEvaluator{}, true);
}

class GasMixture
{
public:
    using Gases = std::vector<std::unique_ptr<Gas>>;
    Gas& ensure_gas(const std::string& name)
    {
        for (const auto& gas : m_gases)
        {
            if (gas->name()==name)
            {
                return *gas;
            }
        }
        return *m_gases.emplace_back(new Gas{name});
    }
    template <class ...Args>
    auto ensure_state(const std::string& gas, Args... args)
    {
        return ensure_gas(gas).ensure_state(args...);
    }
    auto ensure_state(const std::string& gas)
    {
        return &ensure_gas(gas).states();
    }
    void print(std::ostream& os) const
    {
        for (const auto& gas : m_gases)
        {
            gas->printStates(os);
        }
    }
    void checkPopulations() const
    {
        for (const auto& gas : m_gases)
        {
            gas->checkPopulations();
        }
    }
    void evaluateStateDensities() const
    {
        for (const auto& gas : m_gases)
        {
            gas->evaluateStateDensities();
        }
    }
private:
    Gases m_gases;
};

} // namespace experimental
} // namespace loki

#endif // LOKI_CPP_STATECONTAINER_H
