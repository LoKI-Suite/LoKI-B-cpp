/** \file
 *
 *  Functions for calculating enum-values from their string representation.
 *
 *  \author Jan van Dijk
 *  \date   July 2020
 */

#include "LoKI-B/Enumeration.h"
#include "LoKI-B/Log.h"
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

    /** \todo DB: How to get LoKI-B's type? There can be multiple tags, how do these map to LoKI-B's type?
     *  Example: e + N2(X) -> e + N2(A3Su+,v=1-4) has type 'Excitation' in v4, but
     *  type_tags "Electronic", "Vibrational" in v4.
     */
    CollisionType getCollisionTypeFromTypeTagArray(const json_type& type_tags)
    {
        const std::vector<std::string>& tags = type_tags;
        if (tags.size()==0)
        {
            throw std::runtime_error("No type tags found.");
        }
        else if (tags.size()==1)
        {
            if (tags[0]=="Electronic")
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
            if (std::find(tags.begin(),tags.end(),"Ionization")!=tags.end())
            {
                Log<Message>::Notify("Found 'Ionization'. Setting the collision type to Ionization, ignoring other tags.");
                return getCollisionType("Ionization");
            }
            else if (std::find(tags.begin(),tags.end(),"Electronic")!=tags.end())
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
