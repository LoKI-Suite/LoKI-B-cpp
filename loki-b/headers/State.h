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
    // TODO: This struct is temporary!
    using namespace Enumeration;

    struct StateEntry {
        StateType level;
        std::string gasName, e, v, J;
        int16_t charge;
    };

    template <typename TraitType>
    class State {
    protected:
        StateType type;
        std::string e, v, J;
        int16_t charge;

        double energy{0.},
               statisticalWeight{0.},
               population{0.},
               density{0.};

        State(const StateEntry &entry,
                typename Trait<TraitType>::Gas *gas,
                typename Trait<TraitType>::State * parent);

    public:
        typename Trait<TraitType>::State *parent{nullptr};
        typename Trait<TraitType>::Gas *gas{nullptr};

        std::vector<typename Trait<TraitType>::State *> children;

        bool operator>=(const StateEntry &entry);

        void print() const;

        virtual ~State() = default;
    };

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
        std::cout << "charge = " << charge << " e = " << e;
        if (type >= StateType::vibrational) std::cout << " v = " << v;
        if (type >= StateType::rotational) std::cout << " J = " << J;
        std::cout << std::endl;

        for (auto state : children) {
            if (state != nullptr) state->print();
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
