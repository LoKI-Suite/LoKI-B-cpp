//
// Created by daan on 14-5-19.
//

#ifndef LOKI_CPP_GASMIXTURE_H
#define LOKI_CPP_GASMIXTURE_H

// DONE: Fill the Gas, State and Collision structures.
// DONE: Implement the templated approach to gases and states
// DONE: Think of how we can parse the LXCat data. Especially for states, do we need an
//  auxiliary structure like StateEntry or can we immediately create the state?
//  If not then alter the GasMixture class (replace StateEntry by State).
// DONE: Implement LXCat file parsing.
// DONE: Implement Property file / function / value parsing.

#include "LoKI-B/GasMixtureBase.h"
#include "LoKI-B/Collision.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/Traits.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Parse.h"
#include "LoKI-B/json.h"

#include <vector>
#include <set>
#include <fstream>

namespace loki {


    /** The GasMixture class acts as a base class to the EedfGasMixture and future
     *  ChemGasMixture classes. This class is templated to allow the use of trait
     *  classes (classes that define types). Its constructor is protected such
     *  that the class cannot be instantiated as is, although other classes can be
     *  derived from this class. They can either inherit from GasMixture<Boltzmann>
     *  or GasMixture<Chemistry>. The idea is that the first handles EedfGases,
     *  States and Collisions, whereas the latter handles their Chemistry
     *  equivalents.
     */
    template<typename TraitType>
    class GasMixture : public GasMixtureBase
    {
    public:
        using Gas = typename Trait<TraitType>::Gas;
        using State = typename Trait<TraitType>::State;
        using Collision = typename Trait<TraitType>::Collision;
        // Vector of pointers to all the gases in the mixture.
        const std::vector<Gas*>& gases() const { return m_gases; }
    protected:
        /// \todo Still needed?
        Gas* findGas(const std::string& name) { return static_cast<Gas*>(GasMixtureBase::findGas(name)); }
#if 0
        /// \todo Still needed?
        State* findState(const StateEntry& entry) { return static_cast<State*>(GasMixtureBase::findState(entry)); }
#endif
        /** Tries to add a new gas based on a given name and returns a pointer to it. If a
         *  gas with the same name already exists it returns a pointer to that gas.
         */
        Gas* addGas(const std::string &name)
        {
            Gas* gas = findGas(name);
            if (gas)
            {
                return gas;
            }
            else
            {
                gas = new Gas(name);
                GasMixtureBase::addGas(gas);
                m_gases.emplace_back(gas);
                return gas;
            }
        }

        /** This overload accepts only a reference to a StateEntry object. This is the main
         *  function to call when a new state is to be added to the mixture. It will add the
         *  corresponding gas and ancestor states if the do not yet exist, before adding the
         *  state itself. If the state already exists, then a pointer to this state is
         *  returned instead.
         */
        State* addState(const StateEntry &entry)
        {
            Gas* gas = addGas(entry.gasName);
            auto *state = addState(gas, entry);

            for (uint8_t lvl = electronic; lvl < entry.level; ++lvl) {
                state = addState(state, entry);
            }

            return state;
        }
        /** Creates a collision based on a provided CollisionEntry object. The
         *  gases and states involved in the collision are first created and
         *  added to the mixture. Then a Collision object is created and its
         *  pointer is returned.
         */
        // arguments: smth. like "He + e", "->", "He + e", "Elastic"
        Collision* createCollision(const std::string& lhs, const std::string& sep, const std::string& rhs, const std::string& type)
        {

            std::vector <StateEntry> entry_reactants, entry_products;
            std::vector <uint16_t> entry_stoiCoeff;
            Enumeration::CollisionType entry_type;
            bool entry_isReverse;

            if (!Parse::entriesFromString(lhs, entry_reactants))
                Log<LXCatError>::Error(lhs);
            if (!Parse::entriesFromString(rhs, entry_products, &entry_stoiCoeff))
                Log<LXCatError>::Error(rhs);
            entry_type = Enumeration::getCollisionType(type);
            entry_isReverse = (sep[0] == '<');

            std::vector<GasBase::StateBase*> reactants;
            std::vector<GasBase::StateBase*> products;
            std::set<GasBase*> targetGases;

            for (auto &stateEntry : entry_reactants) {
                auto *state = reactants.emplace_back(addState(stateEntry));
                targetGases.insert(&state->gas());
            }

            if (targetGases.size() != 1)
                Log<Message>::Error("Multiple target gases in a single collision.");

            for (auto &stateEntry : entry_products) {
                products.emplace_back(addState(stateEntry));
            }

            return new Collision(entry_type, reactants, products,
                        entry_stoiCoeff, entry_isReverse);
        }
        Collision* createCollision(const json_type& rcnf)
        {

            std::vector <StateEntry> entry_reactants, entry_products;
            std::vector <uint16_t> entry_stoiCoeff;
            Enumeration::CollisionType entry_type;
            bool entry_isReverse;

            if (!Parse::entriesFromJSON(rcnf.at("lhs"), entry_reactants))
                Log<LXCatError>::Error(rcnf.at("lhs").dump(2));
            if (!Parse::entriesFromJSON(rcnf.at("rhs"), entry_products, &entry_stoiCoeff))
                Log<LXCatError>::Error(rcnf.at("rhs").dump(2));
            entry_type = Enumeration::getCollisionType(rcnf.at("type"));

            entry_isReverse = rcnf.at("superelastic");

            std::vector<GasBase::StateBase*> reactants;
            std::vector<GasBase::StateBase*> products;
            std::set<GasBase*> targetGases;

            for (const auto &stateEntry : entry_reactants) {
                auto *state = reactants.emplace_back(addState(stateEntry));
                targetGases.insert(&state->gas());
            }

            if (targetGases.size() != 1)
                Log<Message>::Error("Multiple target gases in a single collision.");

            for (auto &stateEntry : entry_products) {
                products.emplace_back(addState(stateEntry));
            }

            return new Collision(entry_type, reactants, products,
                        entry_stoiCoeff, entry_isReverse);
        }
    private:
        /** This overload accepts a pointer to a gas. It will check if the electronic
         *  ancestor of the state described by the given StateEntry already exists.
         *  If it does, a pointer is returned to this state, and if it does not, the
         *  ancestor is added its pointer returned.
         */

        State* addState(Gas*gas, const StateEntry &entry)
        {
            auto &states = (entry.charge.empty() ? gas->stateTree : gas->ionicStates);
            auto &gas_states_base = (entry.charge.empty() ? static_cast<GasBase*>(gas)->stateBaseTree : static_cast<GasBase*>(gas)->ionicBaseStates);

            auto *state = gas->find(entry);

            if (state == nullptr) {
                state = gas->states.emplace_back(new State(entry, gas)).get();
                gas_states_base.push_back(state);
                /// \todo also add to the mixture's base state vector, once that has added.
                return states.emplace_back(state);
            }

            return state;
        }

        /** This overload accepts a pointer to a state. It will check if one of the
         *  current children of this state is an ancestor of (or equal to) the state
         *  described by the given StateEntry. If such a child exists, its pointer is
         *  returned, if it does not exist then the appropriate state is created and
         *  its pointer returned.
         */
        State* addState(State* parent, const StateEntry &entry)
        {
            auto * state = parent->find(entry);
            if (state == nullptr)
            {
                state = parent->gas()->states.emplace_back(
                        new State(entry, parent->gas(), parent)).get();
                parent->add_child(state);
            }
            return state;
        }
        // Vector of pointers to all the gases in the mixture.
        std::vector<Gas*> m_gases;
    };

} // namespace loki

#endif //LOKI_CPP_GASMIXTURE_H
