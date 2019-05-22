//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_STATE_H
#define LOKI_CPP_STATE_H

#include "Enumeration.h"
#include "Traits.h"

#include <iostream>
#include <string>
#include <vector>

namespace loki {
    using namespace Enumeration;

    // TODO: Move this struct
    struct StateEntry {
        StateType level;
        std::string charge, gasName, e, v, J;

        StateEntry(StateType level, std::string &&gasName, std::string &&charge,
                std::string &&e, std::string &&v, std::string &&J)
                : level(level), charge(std::move(charge)), gasName(std::move(gasName)),
                e(std::move(e)), v(std::move(v)), J(std::move(J)) {}
    };

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

        bool operator>=(const StateEntry &entry);

        void print() const;

        virtual ~State() = default;
    };

    template<typename TraitType>
    State<TraitType>::State(StateType type, typename Trait<TraitType>::Gas *gas,
            std::string e, std::string v, std::string J, std::string charge)
            : type(type), gas(gas), e(std::move(e)), v(std::move(v)), J(std::move(J)), charge(std::move(charge)) {}

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
    void State<TraitType>::print() const {
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
