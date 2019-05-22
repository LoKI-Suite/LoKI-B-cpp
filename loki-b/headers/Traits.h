//
// Created by daan on 17-5-19.
//

#ifndef LOKI_CPP_TRAITS_H
#define LOKI_CPP_TRAITS_H

#include <vector>

namespace loki {
    class Boltzmann;
    class Chemistry;

    class EedfGas;
    class EedfState;
    class EedfCollision;

    class ChemGas;
    class ChemState;
    class ChemCollision;

    template<typename TraitType>
    struct Trait;

// TODO: Add typedefs to structures of the other class (for equivalence pointers).

    template<>
    struct Trait<Boltzmann> {
        typedef EedfGas Gas;
        typedef EedfState State;
        typedef EedfCollision Collision;
        typedef EedfState *Reactants;
    };

    template<>
    struct Trait<Chemistry> {
        typedef ChemGas Gas;
        typedef ChemState State;
        typedef ChemCollision Collision;
        typedef std::vector<ChemState *> Reactants;
    };
}

#endif //LOKI_CPP_TRAITS_H
