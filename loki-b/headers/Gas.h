//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_GAS_H
#define LOKI_CPP_GAS_H

#include "Traits.h"
#include "EedfCollision.h"
#include "Log.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <cmath>

#include "GasBase.h"

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
    class Gas : public GasBase {
    public:
        using State = typename Trait<TraitType>::State;

        // TODO: Do we actually need this first vector?
        // Vector that stores pointers to all states in the system.
        std::vector<std::unique_ptr<State>> states;

        // The stateTree vector stores pointers to the electronic non-ionic states.
        // The ionicStates vector stores pointers to the electronic ionic states.
        std::vector<State*> stateTree, ionicStates;

        /// \todo Needed?
        State* find(const StateEntry& entry) { return static_cast<State*>(GasBase::find(entry)); }
    protected:
        explicit Gas(std::string name) : GasBase(name) {}
    };
}


#endif //LOKI_CPP_GAS_H
