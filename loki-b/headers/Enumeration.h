//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_ENUMERATION_H
#define LOKI_CPP_ENUMERATION_H

#include <cstdint>

namespace loki {
    namespace Enumeration {

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

        enum StatePropertyType : uint8_t {
            direct,
            file,
            function
        };

    }
}

#endif //LOKI_CPP_ENUMERATION_H
