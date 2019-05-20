//
// Created by daan on 17-5-19.
//

#ifndef LOKI_CPP_TRAITS_H
#define LOKI_CPP_TRAITS_H

namespace loki {
    class Boltzmann;
    class Chemistry;

    class EedfGas;
    class EedfState;

    class ChemGas;
    class ChemState;

    template<typename TraitType>
    struct Trait;

// TODO: Add typedefs for Collisions in Trait classes

    template<>
    struct Trait<Boltzmann> {
        typedef EedfGas Gas;
        typedef EedfState State;
    };

    template<>
    struct Trait<Chemistry> {
        typedef ChemGas Gas;
        typedef ChemState State;
    };
}

#endif //LOKI_CPP_TRAITS_H
