//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_ENUMERATION_H
#define LOKI_CPP_ENUMERATION_H

#include <cstdint>
#include <string>
#include <vector>
#include <stdexcept>

namespace loki::Enumeration {

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
        electronic,
        vibrational,
        rotational,
        none
    };

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
    CollisionType getCollisionType(const std::string &str);

}

#endif //LOKI_CPP_ENUMERATION_H
