#ifndef LOKI_CPP_GASBASE_H
#define LOKI_CPP_GASBASE_H

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/StateEntry.h"
#include <iostream>
#include <memory>
#include <vector>

namespace loki
{

class GasBase
{
  public:
    class StateBase
    {
      public:
        using ChildContainer = std::vector<StateBase *>;
        StateBase(const StateEntry &entry, GasBase *gas, StateBase &parent);
        StateBase(GasBase *gas);
        virtual ~StateBase();

        const GasBase &gas() const
        {
            return *m_gas_base;
        }
        GasBase &gas()
        {
            return *m_gas_base;
        }

        const StateBase *parent() const
        {
            return m_parent_base;
        }
        StateBase *parent()
        {
            return m_parent_base;
        }

        const ChildContainer &children() const
        {
            return m_children;
        }
        ChildContainer &children()
        {
            return m_children;
        }

        void printChildren(std::ostream &os) const;
        /* Allows to search the children in order to find a state that is equal
         * to or an ancestor of the state as described by the passed StateEntry
         * object. If this state is present, it will be returned, otherwise a
         * null pointer is returned.
         */
        StateBase *find(const StateEntry &entry);

        /* Returns a vector containing the sibling states of the current state.
         * In case this state is electronic, this vector is obtained via the
         * gas, otherwise it is obtained via the parent.
         */
        const ChildContainer &siblings() const
        {
            if (!parent())
            {
                throw std::runtime_error("Cannot call siblings on the root node.");
            }
            return parent()->m_children;
        }

        /** Verifies that the populations of all child states add up to 1. It also
         *  calls the checkPopulation function on all these states to recursively check
         *  populations.
         */
        void checkPopulations() const;
        /** This operator overload checks whether the current state is equal to
         *  or a parent of a state as described by a given StateEntry object.
         */
        bool operator>=(const StateEntry &entry);
        /* Calculates the density based on the gas fraction, the states population and populations
         * of any parent states. Then it calls evaluateDensity on all child states to recursively
         * evaluate all densities in the state tree.
         */
        void evaluateDensity();
        StateBase *ensure_state(const StateEntry &entry)
        {
            StateBase *state = find(entry);
            if (state)
            {
                return state;
            }
            state = new StateBase(entry, &gas(), *this);
            add_child(state);
            return state;
        }

      protected:
        void add_child(StateBase *child)
        {
            m_children.emplace_back(child);
        }

      private:
        /// \todo Make const?
        GasBase *m_gas_base;
        /// \todo Make const?
        StateBase *m_parent_base;

      public:
        /// \todo make private
        ChildContainer m_children;

      public:
        const StateType type;
        const std::string charge, e, v, J;
        double population;
        double energy;
        double statisticalWeight;
        double density;
    };
    /// \todo Get rid of StateBase, use State exclusively here and elsewhere.
    using State = StateBase;
    explicit GasBase(std::string name);
    virtual ~GasBase();
    /** Prints the (first non-ionic then ionic) electronic states and their
     *  children.
     */
    void print(std::ostream &os) const;
    /* Propagates the gas fractions to the state root and charge levels. If the
     * fraction is >0, the populations of the root node and its neutral child are
     * set to one. Otherwise, populations at these levels remain zero.
     */
    void propagateFraction();
    /* Verifies that the populations of all electronic states adds up to 1. It also
     * calls the checkPopulation function on all these states to recursively check
     * populations.
     */
    void checkPopulations();
    /* Calls evaluateDensity on all electronic states for this gas.
     */
    void evaluateStateDensities();
    /** This member is the State part of the member findState that was previously
     *  completely in GasMixture.h. It differs from member find, which also takes
     *  an Entry and returns a State pointer, but the semantics are not documented
     *  anywhere.
     *  \todo Document the semantics of this function, describe the differences
     *        with find, see if these members can be merged somehow.
     */
    StateBase *findState(const StateEntry &entry);
    State *ensureState(const StateEntry &entry)
    {
        auto *state = get_root().ensure_state(entry);
        // state is now a 'charge state' (container)
        for (uint8_t lvl = StateType::charge; lvl < entry.level; ++lvl)
        {
            state = state->ensure_state(entry);
        }
        return state;
    }

    State &get_root()
    {
        return *m_root;
    }
    const State &get_root() const
    {
        return *m_root;
    }

  public:
    std::unique_ptr<StateBase> m_root;

    const std::string name;
    double mass;
    double harmonicFrequency;
    double anharmonicFrequency;
    double rotationalConstant;
    double electricDipoleMoment;
    double electricQuadrupoleMoment;
    double polarizability;
    double fraction;
};

std::ostream &operator<<(std::ostream &os, const GasBase::StateBase &state);

} // namespace loki

#endif // LOKI_CPP_GASBASE_H
