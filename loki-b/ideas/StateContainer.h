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

#if 0

/// \todo Use this later, instead of the ad-hoc members below (Info base class or member)?
class StateData
{
public:
    double population{0};
    double energy{-1};
    double statisticalWeight{-1};
    double density{0};
};
#endif

class RootInfo
{
public:
    RootInfo(const std::string& gasname)
     : m_gasname(gasname)
    {
    }
    std::string str() const { return m_gasname; }
private:
    const std::string m_gasname;
};

class ChargeInfo
{
public:
    ChargeInfo(const json_type& info) : m_q(info.at("charge").get<std::string>()) {};
    ChargeInfo(const std::string& info) : m_q(info) {};
    bool operator==(const std::string& info) const
    {
        return m_q==info;
    }
    std::string str() const { return "charge="+m_q; }
    double population{0};
private:
    const std::string m_q;
};

class ElectronicInfo
{
public:
    ElectronicInfo(const json_type& info) : m_e(info.at("e").get<std::string>()) {};
    ElectronicInfo(const std::string& info) : m_e(info) {};
    bool operator==(const std::string& info) const
    {
        return m_e==info;
    }
    std::string str() const { return "e="+m_e; }
    double population{0};
private:
    std::string m_e;
};

class VibrationalInfo
{
public:
    VibrationalInfo(const json_type& info) : m_v(info.at("v").get<std::string>()) {};
    VibrationalInfo(const std::string& info) : m_v(info) {};
    bool operator==(const std::string& info) const
    {
        return m_v==info;
    }
    std::string str() const { return "v="+m_v
        + ", g_v=" + std::to_string(statisticalWeight)
        + ", E/eV=" + std::to_string(energy)
        + ", population=" + std::to_string(population); }
    const std::string& v() const { return m_v; }

    // ad hoc public members to test propertyFunction ideas...
    double energy=-1;
    double statisticalWeight=-1;
    double population{0};
private:
    std::string m_v;
};

class RotationalInfo
{
public:
    RotationalInfo(const json_type& info) : m_J(info.at("J").get<std::string>()) {};
    RotationalInfo(const std::string& info) : m_J(info) {};
    bool operator==(const std::string& info) const
    {
        return m_J==info;
    }
    std::string str() const { return "J="+m_J
        + ", g_J=" + std::to_string(statisticalWeight)
        + ", E/eV=" + std::to_string(energy)
        + ", population=" + std::to_string(population); }
    const std::string& J() const { return m_J; }

    // ad hoc public members to test propertyFunction ideas...
    double energy=-1;
    double statisticalWeight=-1;
    double population{0};
private:
    std::string m_J;
};


//enum StateLevel { Charge, Electronic, Vibrational, Rotational };

template <typename ChildT>
class ChildContainer
{
    using ChildrenType = std::vector<std::unique_ptr<ChildT>>;
public:
    using ChildType = ChildT;
    using size_type = typename ChildrenType::size_type;
    using iterator = typename ChildrenType::iterator;
    using const_iterator = typename ChildrenType::const_iterator;

    bool empty() const { return m_children.empty(); }
    size_type size() const { return m_children.size(); }
    iterator begin() { return m_children.begin(); }
    iterator end() { return m_children.end(); }
    const_iterator begin() const { return m_children.begin(); }
    const_iterator end() const { return m_children.end(); }

    void print(std::ostream& os, const std::string& prefix=std::string{}) const
    {
        for (const auto& c : m_children)
        {
            c->print(os,prefix);
        }
    }
    template <class StateT>
    ChildT* add(StateT& parent, const std::string& info)
    {
        return m_children.emplace_back(new ChildType(parent,info)).get();
    }
    ChildType* find(const std::string& info) const
    {
        for (const auto& c : m_children)
        {
            if (*c==info)
            {
                return c.get();
            }
        }
        return nullptr;
    }
    void checkPopulations(const std::string& name) const
    {
        if (m_children.empty())
        {
            return;
        }
        double totalPopulation = 0.;
        for (auto& state : m_children)
        {
            totalPopulation += state->info().population;
            state->checkPopulations();
        }
        if (std::abs(totalPopulation - 1.) > 10. * std::numeric_limits<double>::epsilon())
        {
            throw std::runtime_error("Populations of children of '" + name + "' do not add up to unity.");
        }
    }
private:
    ChildrenType m_children;
};

template <>
class ChildContainer<void>
{
public:
    void print(std::ostream& os, const std::string& prefix) const
    {
    }
    void checkPopulations(const std::string& name) const
    {
    }
};

template <typename StateInfo, typename ParentT, typename StateT, typename ChildT=void>
class State
{
public:
    using Parent = ParentT;
    using ChildContainer = loki::experimental::ChildContainer<ChildT>;

    const Parent& parent() const { return m_parent; }
    const StateInfo& info() const { return m_info; }
    // reconsider:
    StateInfo& info() { return m_info; }
    const ChildContainer& children() const { return m_children; }

    ChildT* ensure_state(const std::string& info)
    {
        ChildT* c = m_children.find(info);
        return c ? c : m_children.add(*static_cast<StateT*>(this),info);
    }
    template <class ...Args>
    auto ensure_state(const std::string& info, Args... args)
    {
        return ensure_state(info)->ensure_state(args...);
    }
    bool operator==(const std::string& info) const
    {
        return m_info==info;
    }
    void print(std::ostream& os, const std::string& prefix=std::string{}) const
    {
        os << prefix << m_info.str() << std::endl;
        m_children.print(os,prefix+' ');
    }
    State(const ParentT& parent, const std::string& name)
    : m_parent(parent),
      m_info(name)
    {
        static_assert(std::is_base_of_v<State,StateT>);
    }
    void checkPopulations() const
    {
        m_children.checkPopulations(m_info.str());
    }
private:
    const Parent& m_parent;
    ChildContainer m_children;
    StateInfo m_info;
};

template <class StateT>
const typename StateT::Parent::ChildContainer& siblings(const StateT& state)
{
    return state.parent().children();
}

class Gas;
class StateRoot;
class StateCharge;
class StateElectronic;
class StateVibrational;
class StateRotational;

class StateRoot : public State<RootInfo,Gas,StateRoot,StateCharge>
{
public:
    StateRoot(Gas& gas);
    const Gas& getGas() const { return parent(); }
};

class StateCharge : public State<ChargeInfo,StateRoot,StateCharge,StateElectronic>
{
public:
    using State::State;
    const Gas& getGas() const { return parent().getGas(); }
};

class StateElectronic : public State<ElectronicInfo,StateCharge,StateElectronic,StateVibrational>
{
public:
    using State::State;
    const Gas& getGas() const { return parent().getGas(); }
};

class StateVibrational : public State<VibrationalInfo,StateElectronic,StateVibrational,StateRotational>
{
public:
    using State::State;
    const Gas& getGas() const { return parent().getGas(); }
};

class StateRotational : public State<RotationalInfo,StateVibrational,StateRotational>
{
public:
    using State::State;
    const Gas& getGas() const { return parent().getGas(); }
};

class Gas
{
public:
    explicit Gas(const std::string& name);
    /** arguments: C [, E [, V [, J ] ] ],
     * where C, E, V and J are info-strings for the charge, electronic,
     * vibrational and rotational aspects of the state. As an example,
     * ensure_state("0","X") creates an eleronic state, while the four-argument
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

    void checkPopulations() const
    {
        std::cout << "Checking populations of gas '" << m_name << std::endl;
        m_states->checkPopulations();
    }

    double harmonicFrequency;
    double rotationalConstant;
private:
    const std::string m_name;
    std::unique_ptr<StateRoot> m_states;
};


// implementation:

StateRoot::StateRoot(Gas& gas)
    : State<RootInfo,Gas,StateRoot,StateCharge>(gas,gas.name())
{
}

inline Gas::Gas(const std::string& name)
    : m_name(name),
    m_states(new StateRoot(*this)), harmonicFrequency(-1), rotationalConstant(-1)
{
}

void Gas::printStates(std::ostream& os) const
{
    m_states->print(os,std::string{});
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
private:
    Gases m_gases;
};

} // namespace experimental
} // namespace loki

#endif // LOKI_CPP_STATECONTAINER_H
