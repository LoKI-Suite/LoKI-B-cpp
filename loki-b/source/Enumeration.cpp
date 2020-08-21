/** \file
 *
 *  Functions for calculating enum-values from their string representation.
 *
 *  \author Jan van Dijk
 *  \date   July 2020
 */

#include "LoKI-B/Enumeration.h"
#include <vector>
#include <utility>
#include <stdexcept>

namespace loki {

namespace {

    /** this template is used for the implementation of the functions
     *  getEedfType, getIonizationOperatorType and getGrowthModelType.
     */
    template <class E>
    E parse_enum_string(const std::string& key, const std::vector<std::pair<std::string,E>>& map)
    {
        for (const auto& e : map)
        {
		if (e.first==key)
                {
                    return e.second;
                }
        }
        throw std::runtime_error("Illegal argument '" + key + "'.");
    }

} // namespace


    EedfType getEedfType(const std::string& str)
    {
        return parse_enum_string<EedfType>(str, {
                { "boltzmann", EedfType::boltzmann },
                { "prescribed", EedfType::prescribed },
                });
    }

    IonizationOperatorType getIonizationOperatorType(const std::string& str)
    {
        return parse_enum_string<IonizationOperatorType>(str, {
                { "conservative", IonizationOperatorType::conservative },
                { "oneTakesAll", IonizationOperatorType::oneTakesAll },
                { "equalSharing", IonizationOperatorType::equalSharing },
                { "usingSDCS", IonizationOperatorType::sdcs }
                });
    }

    GrowthModelType getGrowthModelType(const std::string& str)
    {
        return parse_enum_string<GrowthModelType>(str, {
                { "spatial", GrowthModelType::spatial },
                { "temporal", GrowthModelType::temporal },
                });
    }

    CollisionType getCollisionType(const std::string &str)
    {
        return parse_enum_string<CollisionType>(str, {
                { "Elastic", CollisionType::elastic },
                { "Effective", CollisionType::effective },
                { "Excitation", CollisionType::excitation },
                { "Vibrational", CollisionType::vibrational },
                { "Rotational", CollisionType::rotational },
                { "Ionization", CollisionType::ionization },
                { "Attachment", CollisionType::attachment },
                });
    }

} // namespace loki
