//
// Created by daan on 6-5-19.
//

#ifndef LOKI_CPP_PARSE_H
#define LOKI_CPP_PARSE_H

#include <string>
#include <regex>
#include <iostream>
#include <fstream>

#include "Enumeration.h"
#include "State.h"

namespace loki {
    struct Parse {
        /*
         * The setField function extracts a value from a field in the input file,
         * casts it to the appropriate type and assigns it to the 'value' argument
         * which is passed by reference. It returns a boolean to give an indication
         * whether the operation was successful.
         *
         * Note that this is a template function, and the type of 'value' is arbitrary,
         * a standard string can be converted to most basic types using a stringstream.
         * However, for some basic types (e.g. bool) and custom types, such as enums,
         * we can specialize this function to alter its behaviour.
         *
         * The definitions of these specializations can be found below the Parse struct.
         */

        template<typename T>
        static bool setField(const std::string &sectionContent,
                             const std::string &fieldName, T &value) {

            std::string valueBuffer;

            if (!getFieldValue(sectionContent, fieldName, valueBuffer))
                return false;

            std::stringstream s(valueBuffer);
            s >> value;

            return true;
        }

        /*
         * The getFieldValue function will extract the value of a given field name.
         * From a section in the input file. This function is specifically designed
         * to extract single line values (thus they are not a section or list). The
         * string buffer to hold the value is passed by reference and the function
         * returns a boolean to specify whether the operation was successful.
         *
         * NOTE: This function does not deal with ambiguous field names. E.g. the "isOn"
         * field occurs multiple times in the file. The user is advised to only retrieve
         * fields that are one level above the level of "sectionContent". E.g. retrieve
         * "isOn" when "sectionContent" contains the contents of the "electronKinetics"
         * section.
         */
        static bool getFieldValue(const std::string &sectionContent,
                                  const std::string &fieldName, std::string &valueBuffer) {

            std::regex r(fieldName + R"(:\s*(.*[^\s\n])\s*\n*)");
            std::smatch m;

            if (!std::regex_search(sectionContent, m, r))
                return false;

            valueBuffer = m[1];

            return true;
        }

        /*
         * The getList function retrieves all entries in a list type field, e.g.:
         *
         * dataFiles:
         *   - eedf
         *   - swarmParameters
         *
         * and returns them as a vector of strings (thus {"eedf", "swarmParameters"}
         * in the example). This vector is passed by reference as an argument. and
         * the function returns a boolean to specify whether the operation was
         * successful.
         */
        static bool getList(const std::string &sectionContent,
                            const std::string &fieldName, std::vector<std::string> &container) {

            std::regex r(R"(-\s*(\S+(?: \S+)*)\s*\n*)");

            for (auto it = std::sregex_iterator(sectionContent.begin(), sectionContent.end(), r);
                 it != std::sregex_iterator(); ++it) {

                container.emplace_back(it->str(1));
            }

            return !container.empty();
        }

        /*
         * getSection is a static function that retrieves the contents of a specified section
         * and stores them in the "sectionBuffer" string. Furthermore, it returns a boolean
         * to indicate whether the operation was successful.
         */
        static bool getSection(const std::string &fileContent, const std::string &sectionTitle,
                               std::string &sectionBuffer) {

            // This regular expression finds the level of a specific section. In other
            // words, it finds the number of spaces that precede the section title
            const std::regex reLevel(R"((?:^|\n)( *))" + sectionTitle);
            std::smatch m;

            if (!std::regex_search(fileContent, m, reLevel))
                return false;

            const std::string levelString = m[1];

            // This regular expression matches a specific section in the input file. More accurately,
            // it returns the text in between the specified section and the next section on the same
            // level.
            const std::regex reSection(
                    sectionTitle + R"(:\s*\n*([^]*?)(?:(?:\n+)" + levelString + R"(\w)|(?:\n*\s*$)))");

            if (!std::regex_search(fileContent, m, reSection))
                return false;

            // Restore the original level of the content to allow for further processing.
            sectionBuffer = levelString + "  ";
            sectionBuffer += m[1];

            return true;
        }

        /*
         * removeComments is a static function that takes a string as an argument. It strips the
         * string from comments and returns it.
         */
        static std::string removeComments(const std::string &content) {
            const std::regex reClean(R"(%[^]*?\n)");
            return std::regex_replace(content, reClean, "\n");
        }

        // TODO: comment getFirstValue
        static bool getFirstValue(const std::string &valueString, double &value) {
            if (!getValue(valueString, value))
                return firstValueInRange(valueString, value);

            return true;
        }

        // TODO: comment entriesFromString
        static bool entriesFromString(std::string &statesString, std::vector<StateEntry> &entries,
                                      std::vector<uint16_t> *stoiCoeff = nullptr) {
            std::regex reState(
                    R"((\d*)([A-Z0-9]+)\(([-\+]?)\s*,?\s*([-\+'\[\]/\w]+)\s*(?:,\s*v\s*=\s*([-\+\w]+))?\s*(?:,\s*J\s*=\s*([-\+\d]+))?\s*)");

            std::regex_iterator<std::string::iterator> rit(statesString.begin(),
                                                           statesString.end(),
                                                           reState);
            std::regex_iterator<std::string::iterator> rend;

            if (rit == rend) return false;

            while (rit != rend) {
                Enumeration::StateType stateType;

                if (rit->str(2).empty() || rit->str(4).empty())
                    return false;

                if (rit->str(5).empty()) {
                    stateType = electronic;
                } else if (rit->str(6).empty()) {
                    stateType = vibrational;
                } else {
                    stateType = rotational;
                }

                if (stoiCoeff != nullptr) {
                    if (rit->str(1).empty()) {
                        stoiCoeff->emplace_back(1);
                    } else {
                        std::stringstream ss(rit->str(1));
                        uint16_t coeff;

                        ss >> coeff;
                        stoiCoeff->emplace_back(coeff);
                    }
                }

                entries.emplace_back(stateType, rit->str(2), rit->str(3),
                                     rit->str(4), rit->str(5), rit->str(6));

                ++rit;
            }

            return true;
        }

        // TODO: comment collisionTypeFromString
        static Enumeration::CollisionType collisionTypeFromString(const std::string &collisionTypeString) {
            std::regex rState(R"((Elastic|Effective|Excitation|Vibrational|Rotational|Ionization|Attachment))");
            std::smatch mState;

            if (!std::regex_search(collisionTypeString, mState, rState))
                return Enumeration::CollisionType::none;

            char first = mState.str(0)[0],
                    second = mState.str(0)[1];

            switch (first) {
                case 'E':
                    switch (second) {
                        case 'c':
                            return Enumeration::CollisionType::elastic;
                        case 'e':
                            return Enumeration::CollisionType::effective;
                        default:
                            return Enumeration::CollisionType::excitation;
                    }
                case 'V':
                    return Enumeration::CollisionType::vibrational;
                case 'R':
                    return Enumeration::CollisionType::rotational;
                case 'I':
                    return Enumeration::CollisionType::ionization;
                default:
                    return Enumeration::CollisionType::attachment;
            }
        }

        // TODO: comment rawCrossSectionFromStream
        static void
        rawCrossSectionFromStream(std::vector<std::pair<double, double>> &rawCrossSection, std::ifstream &in) {
            std::string line;

            while(std::getline(in, line)) {
                if (line.substr(0, 2) == "--") break;
            }

            double energy = 0., value = 0.;

            while(std::getline(in, line)) {
                std::stringstream ss(line);
                if (!((ss >> energy) && (ss >> value))) break;\
                rawCrossSection.emplace_back(energy, value);
            }
        }

    private:

        /* -- PARSING RANGES -- */

        // TODO: comment getValue
        static bool getValue(const std::string &valueString, double &value) {
            const std::regex r(R"(\s*(\d+\.?\d*)\s*\n*)");
            std::smatch m;

            if (!std::regex_match(valueString, r)) return false;

            std::stringstream ss(valueString);

            return (bool) (ss >> value);
        }

        // TODO: comment firstValueInRange
        static bool firstValueInRange(const std::string &rangeString, double &value) {
            const std::regex r(R"(\s*(?:(?:logspace\(-?)|(?:linspace\())\s*(\d+\.?\d*)\s*)");
            std::smatch m;

            if (!std::regex_search(rangeString, m, r)) return false;

            std::stringstream ss(m[1]);

            return (bool) (ss >> value);
        }

    };

    /*
     * The setField function needs to behave differently when it is
     * supplied with some specific types. The types for which the
     * function needs to be specialized are:
     *
     * - bool; since these values are supplied through 'true' and
     *   'false' rather than '1' and '0'.
     * - std::string; strictly speaking, this is not necessary,
     *   but in this case we can skip the type cast since the
     *   desired type is already std::string.
     * - std::vector<std::string>; in this case the behaviour
     *   is completely different and the getList function is
     *   used.
     * - enums; whenever a new enum is defined that can occur in
     *   the input file, a new specialization needs to be written.
     *
     * Note that we could have defined a different function for every
     * type (e.g. setBool, setString, ...), however, now we can set
     * any desired value through a single function.
     *
     * Use "inline" with fully specialized templates to refrain from
     * breaking the "One Definition Rule". An alternative is to add
     * a Parse.cpp file and specify the specializations there.
     */
    template<>
    inline
    bool Parse::setField<std::string>(const std::string &sectionContent,
                                      const std::string &fieldName, std::string &value) {

        return getFieldValue(sectionContent, fieldName, value);
    }

    template<>
    inline
    bool Parse::setField<bool>(const std::string &sectionContent,
                               const std::string &fieldName, bool &value) {
        std::string valueBuffer;

        if (!getFieldValue(sectionContent, fieldName, valueBuffer))
            return false;

        std::stringstream s(valueBuffer);
        s >> std::boolalpha >> value;

        return true;
    }

    template<>
    inline
    bool Parse::setField<Enumeration::EedfType>(const std::string &sectionContent,
                                                const std::string &fieldName, Enumeration::EedfType &value) {

        std::string valueBuffer;

        if (!getFieldValue(sectionContent, fieldName, valueBuffer))
            return false;

        if (valueBuffer == "boltzmann") {
            value = Enumeration::EedfType::boltzmann;
        } else if (valueBuffer == "prescribed") {
            value = Enumeration::EedfType::prescribed;
        } else {
            return false;
        }

        return true;
    }

    template<>
    inline
    bool Parse::setField<Enumeration::IonizationOperatorType>(const std::string &sectionContent,
                                                              const std::string &fieldName,
                                                              Enumeration::IonizationOperatorType &value) {

        std::string valueBuffer;

        if (!getFieldValue(sectionContent, fieldName, valueBuffer))
            return false;

        if (valueBuffer == "conservative") {
            value = Enumeration::IonizationOperatorType::conservative;
        } else if (valueBuffer == "oneTakesAll") {
            value = Enumeration::IonizationOperatorType::oneTakesAll;
        } else if (valueBuffer == "equalSharing") {
            value = Enumeration::IonizationOperatorType::equalSharing;
        } else if (valueBuffer == "usingSDCS") {
            value = Enumeration::IonizationOperatorType::sdcs;
        } else {
            return false;
        }

        return true;
    }

    template<>
    inline
    bool Parse::setField<Enumeration::GrowthModelType>(const std::string &sectionContent,
                                                       const std::string &fieldName,
                                                       Enumeration::GrowthModelType &value) {

        std::string valueBuffer;

        if (!getFieldValue(sectionContent, fieldName, valueBuffer))
            return false;

        if (valueBuffer == "spatial") {
            value = Enumeration::GrowthModelType::spatial;
        } else if (valueBuffer == "temporal") {
            value = Enumeration::GrowthModelType::temporal;
        } else {
            return false;
        }

        return true;
    }

    template<>
    inline
    bool Parse::setField<std::vector<std::string>>(const std::string &sectionContent,
                                                   const std::string &fieldName, std::vector<std::string> &value) {

        std::string fieldContent;

        if (!Parse::getSection(sectionContent, fieldName, fieldContent))
            return false;

        return Parse::getList(fieldContent, fieldName, value);

    }
}

#endif //LOKI_CPP_PARSE_H
