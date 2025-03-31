/** \file
 *
 *  Implementation of functions that support enum types.
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

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Log.h"
#include <stdexcept>
#include <utility>
#include <vector>

namespace loki
{

namespace
{

/** this template is used for the implementation of the functions
 *  getEedfType, getIonizationOperatorType and getGrowthModelType.
 */
template <class E>
E parse_enum_string(const std::string &key, const std::vector<std::pair<std::string, E>> &map)
{
    for (const auto &e : map)
    {
        if (e.first == key)
        {
            return e.second;
        }
    }
    throw std::runtime_error("Illegal argument '" + key + "'.");
}

} // namespace

EedfType getEedfType(const std::string &str)
{
    return parse_enum_string<EedfType>(str, {
                                                {"boltzmann", EedfType::boltzmann},
                                                {"prescribed", EedfType::prescribed},
                                            });
}

IonizationOperatorType getIonizationOperatorType(const std::string &str)
{
    return parse_enum_string<IonizationOperatorType>(str, {{"conservative", IonizationOperatorType::conservative},
                                                           {"oneTakesAll", IonizationOperatorType::oneTakesAll},
                                                           {"equalSharing", IonizationOperatorType::equalSharing},
                                                           {"usingSDCS", IonizationOperatorType::sdcs}});
}

GrowthModelType getGrowthModelType(const std::string &str)
{
    return parse_enum_string<GrowthModelType>(str, {
                                                       {"spatial", GrowthModelType::spatial},
                                                       {"temporal", GrowthModelType::temporal},
                                                   });
}

std::string statePropertyName(StatePropertyType type)
{
    switch (type)
    {
        case StatePropertyType::energy: return "energy"; break;
        case StatePropertyType::statisticalWeight: return "statistical weight"; break;
        case StatePropertyType::population: return "population"; break;
    }
    throw std::runtime_error("Illegal StatePropertyType");
}

CollisionType getCollisionType(const std::string &str)
{
    return parse_enum_string<CollisionType>(str, {
                                                     {"Elastic", CollisionType::elastic},
                                                     {"Effective", CollisionType::effective},
                                                     {"Excitation", CollisionType::excitation},
                                                     {"Vibrational", CollisionType::vibrational},
                                                     {"Rotational", CollisionType::rotational},
                                                     {"Ionization", CollisionType::ionization},
                                                     {"Attachment", CollisionType::attachment},
                                                 });
}

/** \todo DB: How to get LoKI-B's type? There can be multiple tags, how do these map to LoKI-B's type?
 *  Example: e + N2(X) -> e + N2(A3Su+,v=1-4) has type 'Excitation' in v4, but
 *  type_tags "Electronic", "Vibrational" in v4.
 */
CollisionType getCollisionTypeFromTypeTagArray(const json_type &type_tags)
{
    const std::vector<std::string> &tags = type_tags;
    if (tags.size() == 0)
    {
        throw std::runtime_error("No type tags found.");
    }
    else if (tags.size() == 1)
    {
        if (tags[0] == "Electronic")
        {
            Log<Message>::Notify("Found 'Electronic'. Setting the collision type to Excitation.");
            return getCollisionType("Excitation");
        }
        else
        {
            return getCollisionType(tags[0]);
        }
    }
    else
    {
        // for now: if Ionization is present, use that.
        //            otherwise: if Electronic is present, use "Excitation" (translate).
        // Otherwise we give up.
        if (std::find(tags.begin(), tags.end(), "Ionization") != tags.end())
        {
            Log<Message>::Notify("Found 'Ionization'. Setting the collision type to Ionization, ignoring other tags.");
            return getCollisionType("Ionization");
        }
        else if (std::find(tags.begin(), tags.end(), "Electronic") != tags.end())
        {
            Log<Message>::Notify("Found 'Electronic'. Setting the collision type to Excitation, ignoring other tags.");
            return getCollisionType("Excitation");
        }
        else
        {
            throw std::runtime_error("type_tags: I do not know how to handle the type_tags array.");
        }
    }
}

} // namespace loki
