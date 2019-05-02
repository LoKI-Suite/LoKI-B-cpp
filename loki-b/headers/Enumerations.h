//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_ENUMERATIONS_H
#define LOKI_CPP_ENUMERATIONS_H

#include <cstdint>

namespace loki {

    enum EedfType : uint8_t {
        boltzmann,
        prescribed
    };

    enum IonizationOperatorType : uint8_t {
        conservative,
        oneTakesAll,
        equalSharing,
        sdcs
    };

    enum GrowthModelType : uint8_t {
        spatial,
        temporal
    };


}

#endif //LOKI_CPP_ENUMERATIONS_H
