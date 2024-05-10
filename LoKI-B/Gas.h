#ifndef LOKI_CPP_GASBASE_H
#define LOKI_CPP_GASBASE_H

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Exports.h"
#include "LoKI-B/GasProperties.h"
#include "LoKI-B/StateEntry.h"
#include <iostream>
#include <memory>
#include <vector>

namespace loki
{

class lokib_export Gas
{
  public:
    class lokib_export State
    {
      public:
        using ChildContainer = std::vector<State *>;
        State(const StateEntry &entry, const Gas *gas, const State &parent);
        State(const Gas *gas);
        State(const State&) = delete;
        State(State&&) = default;
        State& operator=(const State&) = delete;
        State& operator=(State&&) = delete;
	virtual ~State();
        const Gas &gas() const
        {
            return *m_gas;
        }
        const State *parent() const
        {
            return m_parent;
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
        State *find(const StateEntry &entry);

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
        /** Calculates the reduced density delta of this state, then calls this
         *  function on all child states (if any). The reduced density is
         *  obtained by multiplying this state's population with the reduced
         *  density of the parent state. For the root state, the gas fraction
         *  is used as the reduced density of the parent instead.
         *
         *  This function should normally be called only indirectly, by calling
         *  evaluateReducedDensities() on the gas.
         */
        void evaluateReducedDensities();
        State *ensure_state(const StateEntry &entry)
        {
            State *state = find(entry);
            if (state)
            {
                return state;
            }
            state = new State(entry, &gas(), *this);
            m_children.emplace_back(state);
            m_state_deleter.emplace_back(state);
            return state;
        }

        double population() const { return m_population; }
        void setPopulation(double population) { m_population = population; }
        /** The reduced density is the density of this state, divided by the
         *  total particle density of the gases in the mixture. (For a
         *  compound state, the state density is defined as the sum of the
         *  densities of the children.) This quantity is discussed in the
         *  LoKI-B manual \cite Manual_1_0_0 in equations 4a,4b and the
         *  accompanyng text.
         *
         *  \todo It appears that in these definitions, only the neutral gas
         *  considered. The reduced densities of charged states and the electron
         *  are initialized to zero.
         */
        double delta() const { return m_delta; }
      private:
        const Gas *m_gas;
        const State *m_parent;
        ChildContainer m_children;
      public:
        /// \todo Make the following members private, provide accessors if necessary
        const StateType type;
        const std::string charge, e, v, J;
        double energy;
        double statisticalWeight;
      private:
        double m_population;
        double m_delta;
	// to ensure deletion of child states, they are also added to a
	// vector of unique_ptr<const State>.
	std::vector<std::unique_ptr<const State>> m_state_deleter;
    };
    Gas(const GasProperties& gasProps, std::string name);
    virtual ~Gas();
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
    /** Calls evaluateReducedDensities on the root state of this gas. (That
     *  function will recuse into the child states.
     */
    void evaluateReducedDensities();
    /** This member is the State part of the member findState that was previously
     *  completely in GasMixture.h. It differs from member find, which also takes
     *  an Entry and returns a State pointer, but the semantics are not documented
     *  anywhere.
     *  \todo Document the semantics of this function, describe the differences
     *        with find, see if these members can be merged somehow.
     */
    State *findState(const StateEntry &entry);
    State *ensureState(const StateEntry &entry)
    {
        auto *state = get_root().ensure_state(entry);
        // state is now a 'charge state' (container)
        for (uint8_t lvl = StateType::charge; lvl < entry.m_level; ++lvl)
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
    const std::string& name() const { return m_name; }

  private:
    std::unique_ptr<State> m_root;

    const std::string m_name;
  public:
    const double mass;
    const double harmonicFrequency;
    const double anharmonicFrequency;
    const double rotationalConstant;
    //const double electricDipoleMoment;
    const double electricQuadrupoleMoment;
    //const double polarizability;
    const double fraction;
};

lokib_export std::ostream &operator<<(std::ostream &os, const Gas::State &state);

} // namespace loki

#endif // LOKI_CPP_GASBASE_H
