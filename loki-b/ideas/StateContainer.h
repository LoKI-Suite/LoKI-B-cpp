#ifndef LOKI_CPP_STATECONTAINER_H
#define LOKI_CPP_STATECONTAINER_H

#include <vector>
#include <memory>
#include <string>
#include <iostream>

#include "LoKI-B/json.h"

namespace loki {

struct RootInfo
{
    RootInfo(const std::string& dummy)
    {
    }
};

struct ChargeInfo
{
    ChargeInfo(const json_type& info) : m_q(info.at("charge").get<std::string>()) {};
    ChargeInfo(const std::string& info) : m_q(info) {};
    std::string m_q;
};

struct ElectronicInfo
{
    ElectronicInfo(const json_type& info) : m_e(info.at("e").get<std::string>()) {};
    ElectronicInfo(const std::string& info) : m_e(info) {};
    std::string m_e;
};

struct VibrationalInfo
{
    VibrationalInfo(const json_type& info) : m_v(info.at("v").get<std::string>()) {};
    VibrationalInfo(const std::string& info) : m_v(info) {};
    std::string m_v;
};

struct RotationalInfo
{
    RotationalInfo(const json_type& info) : m_J(info.at("J").get<std::string>()) {};
    RotationalInfo(const std::string& info) : m_J(info) {};
    std::string m_J;
};


//enum StateLevel { Charge, Electronic, Vibrational, Rotational };

template <class ParentT>
class ParentSupport
{
public:
    using Parent = ParentT;
    ParentSupport(const Parent& p) : m_parent(p) {}
    const Parent& parent() const { return m_parent; }
private:
    const Parent& m_parent;
};

template <typename ChildT>
class ChildContainer
{
public:
    using ChildType = ChildT;
    using ChildrenType = std::vector<std::unique_ptr<ChildType>>;

    void printChildren(std::ostream& os, const std::string& prefix=std::string{}) const
    {
        if (m_children.empty())
        {
            os << prefix << std::endl;
        }
        else
        {
            for (const auto& c : m_children)
            {
                c->printChildren(os,prefix + " " + c->name());
            }
        }
    }
    template <class StateT>
    ChildT* add_child(StateT& parent, const std::string& info)
    {
        return m_children.emplace_back(new ChildType(parent,info)).get();
    }
    ChildrenType m_children;
};

template <>
class ChildContainer<void>
{
public:
    void printChildren(std::ostream& os, const std::string& prefix) const
    {
        os << prefix << std::endl;
    }
};

template <typename ParentT, typename StateT, typename ChildT=void>
class StateBase : public ParentSupport<ParentT>, public ChildContainer<ChildT>
{
public:
    const std::string& name() const { return m_name; }
    ChildT* add_child(const std::string& info)
    {
        return ChildContainer<ChildT>::add_child(*static_cast<StateT*>(this),info);
    }
    template <class ...Args>
    auto add_child(const std::string& info, Args... args)
    {
        return add_child(info)->add_child(args...);
    }
    friend StateT;
    StateBase(const ParentT& parent, const std::string& name)
    : ParentSupport<ParentT>(parent),
      m_name(name)
    {
    }
    // needed? Rather do this in the State's Info part
    const std::string m_name;
};

class StateRoot;

class Gas;
class StateCharge;

class StateRoot : public StateBase<Gas,StateRoot,StateCharge>
{
public:
    StateRoot(Gas& gas);
};

class StateElectronic;

class StateCharge : public StateBase<StateRoot,StateCharge,StateElectronic>
{
public:
    StateCharge(const Parent& parent, const std::string& q)
    : StateBase<StateRoot,StateCharge,StateElectronic>(parent,"charge="+q),
    m_q{q} {}
private:
    const ChargeInfo m_q;
};

class StateVibrational;

class StateElectronic : public StateBase<StateCharge,StateElectronic,StateVibrational>
{
public:
    StateElectronic(const Parent& parent, const std::string& e)
    : StateBase<StateCharge,StateElectronic,StateVibrational>(parent,"electronic="+e),
    m_e{e} {}
private:
    const ElectronicInfo m_e;
};

class StateRotational;

class StateVibrational : public StateBase<StateElectronic,StateVibrational,StateRotational>
{
public:
    StateVibrational(const Parent& parent, const std::string& v)
    : StateBase<StateElectronic,StateVibrational,StateRotational>(parent,"vibrational="+v),
    m_v{v} {}
private:
    const VibrationalInfo m_v;
};

class StateRotational : public StateBase<StateVibrational,StateRotational>
{
public:
    StateRotational(const Parent& parent, const std::string& J)
    : StateBase<StateVibrational,StateRotational,void>(parent,"rotational="+J),
    m_J{J} {}
private:
    const RotationalInfo m_J;
};

class Gas
{
public:
    explicit Gas(const std::string& name);
    template <class ...Args>
    auto add_state(Args... args)
    {
        return m_states->add_child(args...);
    }
    void printStates(std::ostream& os) const;
    const std::string& name() const { return m_name; }
private:
    const std::string m_name;
    std::unique_ptr<StateRoot> m_states;
};


// implementation:

StateRoot::StateRoot(Gas& gas)
    : StateBase<Gas,StateRoot,StateCharge>(gas,gas.name())
{
}

inline Gas::Gas(const std::string& name)
: m_name(name), m_states(new StateRoot(*this))
{
}

void Gas::printStates(std::ostream& os) const
{
    m_states->printChildren(os,name());
}

} // namespace loki

#endif // LOKI_CPP_STATECONTAINER_H
