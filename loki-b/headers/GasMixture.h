//
// Created by daan on 14-5-19.
//

#ifndef LOKI_CPP_GASMIXTURE_H
#define LOKI_CPP_GASMIXTURE_H

// DONE: Fill the Gas, State and Collision structures.
// DONE: Implement the templated approach to gasses and states
// DONE: Think of how we can parse the LXCat data. Especially for states, do we need an
//  auxiliary structure like StateEntry or can we immediately create the state?
//  If not then alter the GasMixture class (replace StateEntry by State).
// DONE: Implement LXCat file parsing.
// TODO: Implement Property file / function / value parsing.

#include "EedfGas.h"
#include "EedfState.h"
#include "EedfCollision.h"
#include "Collision.h"
#include "Setup.h"
#include "Traits.h"
#include "Log.h"

#include <vector>
#include <set>
#include <fstream>
#include <regex>

namespace loki {
    /* -- GasMixture --
     * The GasMixture class acts as a base class to the EedfGasMixture and future
     * ChemGasMixture classes. This class is templated to allow the use of trait
     * classes (classes that define types). Its constructor is protected such
     * that the class cannot be instantiated as is, although other classes can be
     * derived from this class. They can either inherit from GasMixture<Boltzmann>
     * or GasMixture<Chemistry>. The idea is that the first handles EedfGasses,
     * States and Collisions, whereas the latter handles their Chemistry
     * equivalents.
     */

    template<typename TraitType>
    struct GasMixture {
        // Vector of pointers to all the gasses in the mixture.
        std::vector<typename Trait<TraitType>::Gas *> gasses;

        virtual ~GasMixture() {
            for (auto *gas : gasses) {
                delete gas;
            }
        };

        void print() {
            for (int i = 0; i < gasses.size(); ++i) {
                std::cout << "Gas: " << gasses[i]->name << std::endl;
                gasses[i]->print();
            }
        }

    protected:

        GasMixture() = default;

        /* -- addGas --
         * Tries to add a new gas based on a given name and returns a pointer to it. If a
         * gas with the same name already exists it returns a pointer to that gas.
         */

        typename Trait<TraitType>::Gas *
        addGas(const std::string &name) {
            auto it = std::find_if(gasses.begin(), gasses.end(), [&name](typename Trait<TraitType>::Gas *gas) {
                return (*gas == name);
            });

            if (it == gasses.end()) {
                return gasses.emplace_back(new typename Trait<TraitType>::Gas(name));
            }

            return *it;
        }

        /* -- addState --
         * This overload accepts only a reference to a StateEntry object. This is the main
         * function to call when a new state is to be added to the mixture. It will add the
         * corresponding gas and ancestor states if the do not yet exist, before adding the
         * state itself. If the state already exists, then a pointer to this state is
         * returned instead.
         */

        typename Trait<TraitType>::State *
        addState(const StateEntry &entry) {
            auto *gas = addGas(entry.gasName);

            auto *state = addState(gas, entry);

            for (uint8_t lvl = electronic; lvl < entry.level; ++lvl) {
                state = addState(state, entry);
            }

            return state;
        }

        /* -- addState --
         * This overload accepts a pointer to a gas. It will check if the electronic
         * ancestor of the state described by the given StateEntry already exists.
         * If it does, a pointer is returned to this state, and if it does not, the
         * ancestor is added its pointer returned.
         */

        typename Trait<TraitType>::State *
        addState(typename Trait<TraitType>::Gas *gas, const StateEntry &entry) {
            auto &states = (entry.charge.empty() ? gas->stateTree : gas->ionicStates);

            auto it = std::find_if(states.begin(), states.end(),
                                   [&entry](typename Trait<TraitType>::State *state) {
                                       return *state >= entry;
                                   });

            if (it == states.end()) {
                auto *state = gas->states.emplace_back(new typename Trait<TraitType>::State(entry, gas));
                return states.emplace_back(state);
            }

            return *it;
        }

        /* -- addState --
         * This overload accepts a pointer to a state. It will check if one of the
         * current children of this state is an ancestor of (or equal to) the state
         * described by the given StateEntry. If such a child exists, its pointer is
         * returned, if it does not exist then the appropriate state is created and
         * its pointer returned.
         */

        typename Trait<TraitType>::State *
        addState(typename Trait<TraitType>::State *parent, const StateEntry &entry) {

            auto it = std::find_if(parent->children.begin(), parent->children.end(),
                                   [&entry](typename Trait<TraitType>::State *child) {
                                       return *child >= entry;
                                   });

            if (it == parent->children.end()) {
                auto *state = parent->gas->states.emplace_back(
                        new typename Trait<TraitType>::State(entry, parent->gas, parent));
                return parent->children.emplace_back(state);
            }

            return *it;
        }

        /* -- createCollision --
         * Creates a collision based on a provided CollisionEntry object. The
         * gasses and states involved in the collision are first created and
         * added to the mixture. Then a Collision object is created and its
         * pointer is returned.
         */

        typename Trait<TraitType>::Collision *createCollision(CollisionEntry &entry) {
            std::vector<typename Trait<TraitType>::State *> reactants, products;
            std::set<typename Trait<TraitType>::Gas *> targetGasses;

            for (auto &stateEntry : entry.reactants) {
                auto *state = reactants.emplace_back(addState(stateEntry));
                targetGasses.insert(state->gas);
            }

            if (targetGasses.size() != 1)
                Log<Message>::Error("Multiple target gasses in a single collision.");

            for (auto &stateEntry : entry.products) {
                products.emplace_back(addState(stateEntry));
            }

            return new typename Trait<TraitType>::Collision(entry.type, reactants, products,
                                                            entry.stoiCoeff, entry.isReverse);
        }
    };
}


#endif //LOKI_CPP_GASMIXTURE_H
