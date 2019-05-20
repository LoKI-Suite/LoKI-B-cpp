//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_ENUMERATION_H
#define LOKI_CPP_ENUMERATION_H

#include <cstdint>

namespace loki {
    namespace Enumeration {

        enum class EedfType : uint8_t {
            boltzmann,
            prescribed
        };

        enum class IonizationOperatorType : uint8_t {
            conservative,
            oneTakesAll,
            equalSharing,
            sdcs
        };

        enum class GrowthModelType : uint8_t {
            spatial,
            temporal
        };

        enum class StatePropertyType : uint8_t {
            direct,
            file,
            function
        };

        enum StateType : uint8_t {
            ionic,
            electronic,
            vibrational,
            rotational,
            none
        };

        enum class CollisionType : uint8_t {
            elastic,
            effective,
            excitation,
            vibrational,
            rotational,
            ionization,
            attachment
        };

    }
}

#endif //LOKI_CPP_ENUMERATION_H
