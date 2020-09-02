//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_ENUMERATION_H
#define LOKI_CPP_ENUMERATION_H

#include <cstdint>
#include <string>
#include <iostream>
#include "LoKI-B/json.h"

namespace loki {

    enum class EedfType : uint8_t {
        boltzmann,
        prescribed
    };
    EedfType getEedfType(const std::string& str);

    enum class IonizationOperatorType : uint8_t {
        conservative,
        oneTakesAll,
        equalSharing,
        sdcs
    };
    IonizationOperatorType getIonizationOperatorType(const std::string& str);

    enum class GrowthModelType : uint8_t {
        spatial,
        temporal
    };
    GrowthModelType getGrowthModelType(const std::string& str);

    enum class StatePropertyDataType : uint8_t {
        direct,
        file,
        function
    };

    enum class StatePropertyType : uint8_t {
        energy,
        statisticalWeight,
        population
    };

    enum StateType : uint8_t {
        root,
        charge,
        electronic,
        vibrational,
        rotational,
        none
    };
    inline std::ostream& operator<<(std::ostream& os, StateType type)
    {
        switch(type)
        {
        case root: os << "root"; break;
        case charge: os << "charge"; break;
        case electronic: os << "electronic"; break;
        case vibrational: os << "vibrational"; break;
        case rotational: os << "rotational"; break;
        case none: os << "none"; break;
        }
        return os;
    }

    enum class CollisionType : uint8_t {
        effective,
        elastic,
        excitation,
        vibrational,
        rotational,
        ionization,
        attachment,
        size,
        none
    };
    CollisionType getCollisionType(const std::string& str);
    CollisionType getCollisionTypeFromTypeTagArray(const json_type& type_tags);

} // namespace loki

#endif //LOKI_CPP_ENUMERATION_H
