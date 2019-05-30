//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_GAS_H
#define LOKI_CPP_GAS_H

#include "Traits.h"
#include "EedfCollision.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace loki {
    /* -- Gas --
     * This class acts as a base class for both the EedfGas and the future ChemGas. It
     * is templated to allow the use of trait classes (classes that define types). The
     * constructor of this class is protected such that the user cannot instantiate it
     * as is. However, other classes can inherit from this class. In this case we can
     * either inherit Gas<Boltzmann> or Gas<Chemistry>. The difference here is that
     * the Gas<Boltzmann> class has vectors to pointers of EedfStates whereas
     * Gas<Chemistry> has vectors to pointers of ChemStates.
     *
     * The reason that this class is introduced in this format is that an EedfGas and
     * a ChemGas share a lot of properties. On top of that they share the same structure
     * i.e. they both contain vectors to states, however the type of the states differs
     * between the two, which is fixed by introducing the trait classes. In this way
     * we can implement some of the functionality, such as adding states to a gas,
     * only once instead of separately for both classes.
     */

    template <typename TraitType>
    class Gas {
    public:
        const std::string name;
        double mass{-1},
               harmonicFrequency{-1},
               anharmonicFrequency{-1},
               rotationalConstant{-1},
               electricDipoleMoment{-1},
               electricQuadrupoleMoment{-1},
               polarizability{-1},
               fraction{-1};

        // TODO: Do we actually need this first vector?
        // Vector that stores pointers to all states in the system.
        std::vector<typename Trait<TraitType>::State *> states;

        // The stateTree vector stores pointers to the electronic non-ionic states.
        // The ionicStates vector stores pointers to the electronic ionic states.
        std::vector<typename Trait<TraitType>::State *> stateTree, ionicStates;

        // TODO: comment evaluateStateDensities

        void evaluateStateDensities() {
            for (auto *state : stateTree)
                state->evaluateDensity();

            for (auto *state : ionicStates)
                state->evaluateDensity();
        }

        // TODO: comment find

        typename Trait<TraitType>::State *
        find(const StateEntry &entry) {
            auto &childStates = (entry.charge.empty() ? stateTree : ionicStates);

            auto it = std::find_if(childStates.begin(), childStates.end(),
                                   [&entry](typename Trait<TraitType>::State *state) {
                                       return *state >= entry;
                                   });

            if (it == childStates.end()) {
                return nullptr;
            }

            return *it;
        }

        void checkPopulations() {
            double totalPopulation = 0.;

            for (auto *state : stateTree) {
                totalPopulation += state->population;
                state->evaluatePopulations();
            }
            for (auto *state : ionicStates) {
                totalPopulation += state->population;
                state->evaluatePopulations();
            }

            if (std::abs(totalPopulation - 1.) > 10. * std::numeric_limits<double>::epsilon())
                Log<ChildrenPopulationError>::Error(name);
        }

    protected:
        explicit Gas(std::string name) : name(std::move(name)) {}

    public:
        /* -- == --
         * This overload allows a Gas object to be compared to a string. If the
         * string is equal to the name of the gas this function will return true.
         */

        bool operator==(const std::string &otherName) {
            return name == otherName;
        }

        /* -- print --
         * Prints the (first non-ionic then ionic) electronic states and their
         * children.
         */

        void print() const {
            for (uint16_t i = 0; i < stateTree.size(); ++i) {
                std::cout << *stateTree[i] << std::endl;
                stateTree[i]->printChildren();
            }
            for (uint16_t i = 0; i < ionicStates.size(); ++i) {
                std::cout << *ionicStates[i] << std::endl;
                ionicStates[i]->printChildren();
            }
        }

        ~Gas() {
            for (auto *state : states)
                delete state;
        }
    };
}


#endif //LOKI_CPP_GAS_H
