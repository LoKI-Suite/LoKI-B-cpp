#ifndef LOKI_CPP_STATECONTAINER_H
#define LOKI_CPP_STATECONTAINER_H

#include "ideas/Container.h"
#include <limits>

namespace loki {
namespace experimental {

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
// Root has the Gas as parent (todo: reconsider that), StateCharge as child type.

class StateRoot : public Node<RootInfo,Gas,StateRoot,StateCharge>
{
public:
    StateRoot(Gas& gas);
    const Gas& getGas() const { return parent(); }
    /** \todo see if we can provide level() in a better way without having this member in every State type below
     *  Same for getGas(). Make this a free function template?
     */
    static constexpr unsigned level() { return 0; }
};

class StateCharge : public Node<ChargeInfo,StateRoot,StateCharge,StateElectronic>
{
public:
    using Node::Node;
    const Gas& getGas() const { return parent().getGas(); }
    static constexpr unsigned level() { return Parent::level()+1; }
};

class StateElectronic : public Node<ElectronicInfo,StateCharge,StateElectronic,StateVibrational>
{
public:
    using Node::Node;
    const Gas& getGas() const { return parent().getGas(); }
    static constexpr unsigned level() { return Parent::level()+1; }
};

class StateVibrational : public Node<VibrationalInfo,StateElectronic,StateVibrational,StateRotational>
{
public:
    using Node::Node;
    const Gas& getGas() const { return parent().getGas(); }
    static constexpr unsigned level() { return Parent::level()+1; }
};

class StateRotational : public Node<RotationalInfo,StateVibrational,StateRotational,void>
{
public:
    using Node::Node;
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
    : Node<RootInfo,Gas,StateRoot,StateCharge>(gas,gas.name())
{
}

inline Gas::Gas(const std::string& name)
    :
    harmonicFrequency(-1),
    rotationalConstant(-1),
    fraction(1),
    m_name(name),
    m_states(new StateRoot(*this))
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
