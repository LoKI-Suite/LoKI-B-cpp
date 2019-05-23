//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_STATE_H
#define LOKI_CPP_STATE_H

#include "InputStructures.h"
#include "Enumeration.h"
#include "Traits.h"

#include <iostream>
#include <string>
#include <vector>

namespace loki {
    using namespace Enumeration;

    /* -- State --
     * State is a base class to the EedfState and ChemState classes. It is templated
     * to allow for the use of trait classes (classes that define types). The
     * constructor of this class is protected such that the user cannot instantiate
     * this class as is. However, other classes can be derived from this class. They
     * can either inherit from State<Boltzmann> or State<Chemistry>. The difference
     * here is that a State<Boltzmann> has EedfStates as its parent and children,
     * and a pointer to and EedfGas, whereas a State<Chemistry> has ChemStates as
     * its parent and children and a pointer to a ChemGas.
     *
     * The idea behind the implementation of the State base class is similar to that
     * of the Gas base class.
     */

    template <typename TraitType>
    class State {
    protected:
        StateType type;
        std::string e, v, J, charge;

        double energy{0.},
               statisticalWeight{0.},
               population{0.},
               density{0.};

        State(const StateEntry &entry,
                typename Trait<TraitType>::Gas *gas,
                typename Trait<TraitType>::State * parent);
        State(StateType type, typename Trait<TraitType>::Gas *gas,
                std::string e, std::string v, std::string J, std::string charge);

    public:
        typename Trait<TraitType>::State *parent{nullptr};
        typename Trait<TraitType>::Gas *gas{nullptr};

        std::vector<typename Trait<TraitType>::State *> children;

        /* -- >= --
         * This operator overload checks whether the current state is equal to
         * or a parent of a state as described by a given StateEntry object.
         */

        bool operator>=(const StateEntry &entry);

        void printChildren() const;

        virtual ~State() = default;
    };

    template<typename TraitType>
    State<TraitType>::State(StateType type, typename Trait<TraitType>::Gas *gas,
            std::string e, std::string v, std::string J, std::string charge)
            : type(type), gas(gas), e(std::move(e)), v(std::move(v)), J(std::move(J)), charge(std::move(charge)) {}

    /* -- Constructor --
     * Instantiates a State object as presented by a StateEntry object. Take note that
     * the type of the state is inferred from its parent and not from the StateEntry.
     * This is done to allow one StateEntry object to be used to create a state and its
     * ancestors if they do not yet exist.
     */

    template<typename TraitType>
    State<TraitType>::State(const StateEntry &entry,
            typename Trait<TraitType>::Gas *gas, typename Trait<TraitType>::State *parent)
            : e(entry.e), v(entry.v), J(entry.J), charge(entry.charge), parent(parent),
              gas(gas) {

        if (parent == nullptr)
            type = electronic;
        else
            type = (StateType)(parent->type + 1);

    }

    template<typename TraitType>
    void State<TraitType>::printChildren() const {
        std::string space = "  ";

        for (uint8_t i = electronic; i < type; ++i) space.append("  ");

        for (auto state : children) {
            if (state != nullptr) {
                std::cout << space << *state << std::endl;
                state->print();
            }
        }
    }

    template<typename TraitType>
    bool State<TraitType>::operator>=(const StateEntry &entry) {
        {
            if (type > entry.level || charge != entry.charge) return false;

            if (e == entry.e) {
                if (type == electronic) return true;
            } else {
                return false;
            }

            if (v == entry.v) {
                if (type == vibrational) return true;
            } else {
                return false;
            }

            if (J == entry.J) {
                if (type == rotational) return true;
            } else {
                return false;
            }

            return true;
        }
    }
}

#endif //LOKI_CPP_STATE_H
