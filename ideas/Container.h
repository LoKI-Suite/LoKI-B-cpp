#ifndef LOKI_CPP_CONTAINER_H
#define LOKI_CPP_CONTAINER_H

#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <type_traits>

namespace loki {
namespace experimental {

template <class ChildT>
using ChildContainer = std::vector<std::unique_ptr<ChildT>>;

/* Template arguments:
 *
 *  - Node data (other than parent, children etc.)
 *  - The Parent, Current and Child types. NodeT must be derived from
 *    Node<NodeInfoT,ParentT,NodeT,ChildT>. This is checked with a
 *    static_assert in the constructor.
 *
 *  NOTE the specialization for ChildT=void, which terminates the node sequence.
 */
template <typename NodeInfoT, typename ParentT, typename NodeT, typename ChildT>
class Node : public NodeInfoT
{
public:
    using NodeInfo = NodeInfoT;
    using Parent = ParentT;
    using Child = ChildT;
    using ChildContainer = loki::experimental::ChildContainer<Child>;

    using size_type = typename ChildContainer::size_type;
    using iterator = typename ChildContainer::iterator;
    using const_iterator = typename ChildContainer::const_iterator;

    Node(const Parent& parent, const std::string& name)
    : NodeInfo(name), m_parent(parent)
    {
        static_assert(std::is_base_of_v<Node,NodeT>);
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
        // create child node 'info' (when necessary),
        // call ensure_state on that object with the remaining args.
        return ensure_state(info)->ensure_state(args...);
    }
    Child* find(const std::string& info) const
    {
        for (const auto& c : m_children)
        {
            /// \todo This will call NodeInfo::operator==(info), which is a bit too subtle.
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
        return m_children.emplace_back(new Child(*static_cast<const NodeT*>(this),info)).get();
    }
    const Parent& m_parent;
    ChildContainer m_children;
};

/* Specialization for ChildT=void. This does not have/support child nodes.
 */
template <typename NodeInfoT, typename ParentT, typename NodeT>
class Node<NodeInfoT,ParentT,NodeT,void> : public NodeInfoT
{
public:
    using NodeInfo = NodeInfoT;
    using Parent = ParentT;
    using Child = void;

    Node(const Parent& parent, const std::string& name)
    : NodeInfo(name), m_parent(parent)
    {
        static_assert(std::is_base_of_v<Node,NodeT>);
    }
    const Parent& parent() const { return m_parent; }
private:
    const Parent& m_parent;
};

/** \todo See if we can be smarter about the various overloads, there
 *  are a lot if also Operator can be a const/non-const reference.
 *  Pass all arguments by value and let the caller use std::ref, std::cref?
 */

template <class NodeT, class Operation, std::enable_if_t<!std::is_same_v<typename NodeT::Child,void>>* =nullptr>
void apply(const NodeT& node, const Operation& operation, bool recurse)
{
    operation(node);
    if (recurse)
    {
        for (const auto& c : node)
        {
            apply(*c,operation,recurse);
        }
    }
}

// version for Child==void
template <class NodeT, class Operation, std::enable_if_t<std::is_same_v<typename NodeT::Child,void>>* =nullptr>
void apply(const NodeT& node, const Operation& operation, bool recurse)
{
    operation(node);
}

template <class NodeT, class Operation, std::enable_if_t<!std::is_same_v<typename NodeT::Child,void>>* =nullptr>
void apply(NodeT& node, const Operation& operation, bool recurse)
{
    if (recurse)
    {
        for (const auto& c : node)
        {
            apply(*c,operation,recurse);
        }
    }
}

// version for Child==void
template <class NodeT, class Operation, std::enable_if_t<std::is_same_v<typename NodeT::Child,void>>* =nullptr>
void apply(NodeT& node, const Operation& operation, bool recurse)
{
    operation(node);
}

} // namespace experimental
} // namespace loki

#endif // LOKI_CPP_CONTAINER_H

