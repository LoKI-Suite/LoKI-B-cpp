/** \file
 *
 *  Declarations of various enum types and associated functions.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   2 May 2019
 */

#ifndef LOKI_CPP_ENUMERATION_H
#define LOKI_CPP_ENUMERATION_H

#include "LoKI-B/json.h"
#include <cstdint>
#include <iostream>
#include <string>

namespace loki
{

enum class EedfType : uint8_t
{
    boltzmann,
    prescribed
};
EedfType getEedfType(const std::string &str);

enum class IonizationOperatorType : uint8_t
{
    conservative,
    oneTakesAll,
    equalSharing,
    sdcs
};
IonizationOperatorType getIonizationOperatorType(const std::string &str);

enum class GrowthModelType : uint8_t
{
    spatial,
    temporal
};
GrowthModelType getGrowthModelType(const std::string &str);

enum class StatePropertyType : uint8_t
{
    energy,
    statisticalWeight,
    population
};
std::string statePropertyName(StatePropertyType type);

enum StateType : uint8_t
{
    root,
    charge,
    electronic,
    vibrational,
    rotational,
    none
};
inline std::ostream &operator<<(std::ostream &os, StateType type)
{
    switch (type)
    {
    case root:
        os << "root";
        break;
    case charge:
        os << "charge";
        break;
    case electronic:
        os << "electronic";
        break;
    case vibrational:
        os << "vibrational";
        break;
    case rotational:
        os << "rotational";
        break;
    case none:
        os << "none";
        break;
    }
    return os;
}

enum class CollisionType : uint8_t
{
    effective,
    elastic,
    excitation,
    vibrational,
    rotational,
    ionization,
    attachment,
    size
};
CollisionType getCollisionType(const std::string &str);
CollisionType getCollisionTypeFromTypeTagArray(const json_type &type_tags);

#ifdef __cpp_lib_to_underlying

// make std::to_underlying (from <utility> available in the loki namespace
using std::to_underlying;

#else

/** Provide to_underlying for enum-values in case std::to_underlying is not
 *  available. See https://en.cppreference.com/w/cpp/utility/to_underlying
 *  for std::to_underlying, https://stackoverflow.com/questions/8357240
 *  (the comment trail) for this alternative implementation.
 */
template <typename E>
constexpr auto to_underlying(E e) noexcept
{
    return static_cast<std::underlying_type_t<E>>(e);
}

#endif

} // namespace loki

#endif // LOKI_CPP_ENUMERATION_H
