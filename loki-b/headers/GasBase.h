#ifndef LOKI_CPP_GASBASE_H
#define LOKI_CPP_GASBASE_H

#include "InputStructures.h"
#include "Enumeration.h"

namespace loki {

class GasBase
{
public:
    class StateBase
    {
    public:
        using StateType = Enumeration::StateType;
        StateBase(const StateEntry &entry, StateType type, GasBase* gas, StateBase* parent);
        StateBase(StateType type, std::string e, std::string v, std::string J, std::string charge, GasBase* gas);
        virtual ~StateBase();
        const StateType type;
        const std::string e, v, J, charge;
        void printChildren(std::ostream& os) const;
        /* Allows to search the children in order to find a state that is equal
         * to or an ancestor of the state as described by the passed StateEntry
         * object. If this state is present, it will be returned, otherwise a
         * null pointer is returned.
         */
        StateBase* find(const StateEntry &entry);
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

        const GasBase& gas_base() const { return *m_gas_base; }
        GasBase& gas_base() { return *m_gas_base; }
        const StateBase* parent_base() const { return m_parent_base; }
        StateBase* parent_base() { return m_parent_base; }
        /* Returns a vector containing the sibling states of the current state.
         * In case this state is electronic, this vector is obtained via the
         * gas, otherwise it is obtained via the parent.
         */
        std::vector<StateBase*> &siblings_base()
        {
            if (type == electronic)
            {
                return charge.empty() ? gas_base().stateBaseTree : gas_base().ionicBaseStates;
            }
            /// \todo assert/test the parent exists
            return parent_base()->m_children;
        }
    protected:
        void add_child(StateBase* child) { m_children.push_back(child); }
    private:
        /// \todo Make const?
        GasBase* m_gas_base;
        /// \todo Make const?
        StateBase* m_parent_base;
    public:
        /// \todo make private
        std::vector<StateBase *> m_children;
    public:
        double population;
        double energy;
        double statisticalWeight;
        double density;
    };
    explicit GasBase(std::string name);
    virtual ~GasBase(){}
    /** Prints the (first non-ionic then ionic) electronic states and their
     *  children.
     */
    void print(std::ostream& os) const;
    /* Verifies that the populations of all electronic states adds up to 1. It also
     * calls the checkPopulation function on all these states to recursively check
     * populations.
     */
    void checkPopulations();
    /* Calls evaluateDensity on all electronic states for this gas.
     */
    void evaluateStateDensities();
    /* Allows to search the electronic states in order to find a state that
     * is equal to or an ancestor of the state as described by the passed
     * StateEntry object. If this state is present, it will be returned,
     * otherwise a null pointer is returned.
     */
    StateBase* find(const StateEntry& entry);
    /** This member is the State part of the member findState that was previously
     *  completely in GasMixture.h. It differs from member find, which also takes
     *  an Entry and returns a State pointer, but the semantics are not documented
     *  anywhere.
     *  \todo Document the semantics of this function, describe the differences
     *        with find, see if these members can be merged somehow.
     */
    StateBase* findState(const StateEntry& entry);

    // The stateTree vector stores pointers to the electronic non-ionic states.
    // The ionicStates vector stores pointers to the electronic ionic states.
    std::vector<StateBase*> stateBaseTree, ionicBaseStates;

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

