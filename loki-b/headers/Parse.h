//
// Created by daan on 6-5-19.
//

#ifndef LOKI_CPP_PARSE_H
#define LOKI_CPP_PARSE_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <set>

#include "Enumeration.h"
#include "InputStructures.h"
#include "JobSystem.h"
#include "StandardPaths.h"
#include "json.h"

namespace loki
{
struct Parse
{
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

    template <typename T>
    static bool setField(const std::string &sectionContent, const std::string &fieldName, T &value)
    {

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
    static bool getFieldValue(const std::string &sectionContent, const std::string &fieldName, std::string &valueBuffer)
    {
        const std::regex r(fieldName + R"(:\s*(.*[^\s\n])\s*\n*)");
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
    static bool getList(const std::string &sectionContent, const std::string &fieldName,
                        std::vector<std::string> &container)
    {

        static const std::regex r(R"(-\s*(\S+(?: \S+)*)\s*\n*)");

        for (auto it = std::sregex_iterator(sectionContent.begin(), sectionContent.end(), r);
             it != std::sregex_iterator(); ++it)
        {

            container.emplace_back(it->str(1));
        }

        return !container.empty();
    }

    /*
     * getSection is a static function that retrieves the contents of a specified section
     * and stores them in the "sectionBuffer" string. Furthermore, it returns a boolean
     * to indicate whether the operation was successful.
     */
    static bool getSection(const std::string &fileContent, const std::string &sectionTitle, std::string &sectionBuffer)
    {

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
#ifdef _MSVC
        // MSVC: $ matches before \n and at the end of input
        const std::regex reSection(sectionTitle + R"(:\s*\n*([^]*?)(?:(?:\n+)" + levelString + R"(\w)|(?:$\n$)))");
#else
        // OTHER: $ matches only end of input
        const std::regex reSection(sectionTitle + R"(:\s*\n*([^]*?)(?:(?:\n+)" + levelString + R"(\w)|(?:\n*$)))");
#endif

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
    static std::string removeComments(const std::string &content)
    {
        std::string content_clean;

        static const std::regex reLine(R"(\n\s*%[^\n]*)");
        content_clean = std::regex_replace(content, reLine, "");

        static const std::regex reClean(R"(%[^]*?(?:\n|$))");
        return std::regex_replace(content_clean, reClean, "\n");
    }

    /*  SECOND PARSING PHASE  */

    /* -- getFirstValue --
     * Takes a string that can either be a range (linspace/logspace) or a scalar, and
     * parses it into a double. When the string concerns a range, the first value in
     * the range is returned.
     */

    static bool getFirstValue(const std::string &valueString, double &value)
    {
        if (!getValue(valueString, value))
            return firstValueInRange(valueString, value);

        return true;
    }
    static bool getFirstValue(const json_type &entry, double &value)
    {
        if (entry.type()==json_type::value_t::string)
        {
            value = firstValueInRange(entry, value);
        }
        else
        {
            value = entry.get<double>();
        }
        return true;
    }

    /* -- entriesFromString --
     * Accepts a string containing the LHS or RHS of a collision equation. The states
     * in this expression are then parsed into a vector of StateEntry objects, which
     * is passed by reference. Furthermore, the user can supply a pointer to a vector
     * in which the stoichiometric coefficients of the states in this collision are
     * then stored. If a null pointer is passed, these coefficients are not stored.
     */

    static StateEntry entryFromJSON(const json_type& cnf)
    {
            const std::string gasName = cnf.at("particle");
            const int charge_int = cnf.at("charge").get<int>();
            const std::string charge = charge_int ? std::to_string(charge_int) : std::string{};
            const json_type& descr = cnf.at("descriptor");
            if (descr.contains("states"))
            {
                // We have an array of state objects, instead of a single one.
                // loki-b expects a single string for "e", "v" and "J" (the latter
                // can be empty. We need to 'pessimize' the input by canoncatenating
                // the info into single string. The question here is what loki-b
                // supports. Is it allowed, for example, that the "e" fields are
                // different? For now we assume that everything is possible and we
                // simply concatenate all unique "e"'s while for "v" and "J" we assume
                // that the entries form a continuous value-range.
                std::set<std::string> e_vals;
                std::set<unsigned> v_vals;
                std::set<unsigned> J_vals;
                for (json_type::const_iterator s = descr.at("states").begin(); s!= descr.at("states").end(); ++s)
                {
                    e_vals.insert(s->at("e").get<std::string>());
                    if (s->contains("v"))
                    {
                        v_vals.insert(s->at("v").get<int>());
                    }
                    if (s->contains("J"))
                    {
                        J_vals.insert(s->at("J").get<int>());
                    }
                }
                std::string e;
                if (e_vals.size()!=1)
                {
                    throw std::runtime_error("Expected a unique electronic state identifier.");
                }
                else
                {
                    e = *e_vals.begin();
                }
                std::string v;
                if (v_vals.size()==1)
                {
                    v = *v_vals.begin();
                }
                else if (v_vals.size()>1)
                {
                    int nv = *v_vals.rbegin()+1-*v_vals.begin();
                    if (nv!=v_vals.size())
                    {
                        throw std::runtime_error("Expected a contiguous v-range.");
                    }
                    v = std::to_string(*v_vals.begin()) + '-' + std::to_string(*v_vals.rbegin());
                }
                std::string J;
                if (J_vals.size()==1)
                {
                    J = *J_vals.begin();
                }
                else if (J_vals.size()>1)
                {
                    int nJ = *J_vals.rbegin()+1-*J_vals.begin();
                    if (nJ!=J_vals.size())
                    {
                        throw std::runtime_error("Expected a contiguous J-range.");
                    }
                    J = std::to_string(*J_vals.begin()) + '-' + std::to_string(*J_vals.rbegin());
                }
                Enumeration::StateType stateType
                    = J.empty()==false ? rotational
                    : v.empty()==false ? vibrational
                    : electronic;
                return StateEntry{stateType,gasName,charge,e,v,J};
            }
            else
            {
                const std::string e{descr.at("e").get<std::string>()};
                const std::string v{descr.contains("v") ? (
                    descr.at("v").type()==json_type::value_t::string
                        ? descr.at("v").get<std::string>()
                        : std::to_string(descr.at("v").get<int>())
                ) : std::string{} };
                const std::string J{descr.contains("J") ? std::to_string(descr.at("J").get<int>()) : std::string{} };
                /** \todo Check the precise semantics of the next line. Is it possible that J
                 *        is specified, but not v? Is e always specified?
                 */
                Enumeration::StateType stateType
                    = descr.contains("J") ? rotational
                    : descr.contains("v") ? vibrational
                    : electronic;
                return StateEntry{stateType,gasName,charge,e,v,J};
            }
    }
    static bool entriesFromJSON(const json_type& cnf, std::vector<StateEntry> &entries,
                                  std::vector<uint16_t> *stoiCoeff = nullptr)
    {
        for (json_type::const_iterator i=cnf.begin(); i!=cnf.end(); ++i)
        {
            if (i->at("particle").get<std::string>() == "e")
            {
                continue;
            }
            /** \todo It appeaes that the present JSON simply repeats the particle
             *        object when it appears more than once. Then the stoichiometric
             *        coefficient of each entry will be one, and we hope that 'e + e'
             *        will be handled the same way as '2 e' (JvD).
             */
            entries.push_back(entryFromJSON(*i));
            if (stoiCoeff)
            {
                    stoiCoeff->push_back(1);
            }
        }
        return true;
    }
    static bool entriesFromString(const std::string &statesString, std::vector<StateEntry> &entries,
                                  std::vector<uint16_t> *stoiCoeff = nullptr)
    {
        static const std::regex reState(
            R"((\d*)([A-Za-z][A-Za-z0-9]*)\(([-\+]?)\s*,?\s*([-\+'\[\]/\w]+)\s*(?:,\s*v\s*=\s*([-\+\w]+))?\s*(?:,\s*J\s*=\s*([-\+\d]+))?\s*)");

        std::regex_iterator<std::string::const_iterator> rit(statesString.begin(), statesString.end(), reState);
        std::regex_iterator<std::string::const_iterator> rend;

        if (rit == rend)
            return false;

        while (rit != rend)
        {
            Enumeration::StateType stateType;

            if (rit->str(2).empty() || rit->str(4).empty())
                return false;

            if (rit->str(5).empty())
            {
                stateType = electronic;
            }
            else if (rit->str(6).empty())
            {
                stateType = vibrational;
            }
            else
            {
                stateType = rotational;
            }

            if (stoiCoeff != nullptr)
            {
                if (rit->str(1).empty())
                {
                    stoiCoeff->emplace_back(1);
                }
                else
                {
                    std::stringstream ss(rit->str(1));
                    uint16_t coeff;

                    ss >> coeff;
                    stoiCoeff->emplace_back(coeff);
                }
            }

            entries.emplace_back(stateType, rit->str(2), rit->str(3), rit->str(4), rit->str(5), rit->str(6));
            ++rit;
        }

        return true;
    }

    /* -- propertyStateFromString --
     * Extracts a StateEntry object from a given string and returns it. Note that this function
     * is specifically used when loading state properties, since then the states can contain
     * wild card characters.
     */

    static StateEntry propertyStateFromString(const std::string &propertyString)
    {
        static const std::regex reState(
            R"(([A-Za-z][A-Za-z0-9]*)\(([-\+]?)\s*,?\s*([-\+'\[\]/\w\*]+)\s*(?:,\s*v\s*=\s*([-\+\w\*]+))?\s*(?:,\s*J\s*=\s*([-\+\d\*]+))?\s*)");
        std::smatch m;

        if (!std::regex_search(propertyString, m, reState))
            return {};

        if (m.str(1).empty() || m.str(3).empty())
            return {};

        Enumeration::StateType stateType;

        if (m.str(4).empty())
        {
            stateType = electronic;
        }
        else if (m.str(5).empty())
        {
            stateType = vibrational;
        }
        else
        {
            stateType = rotational;
        }
        return {stateType, m.str(1), m.str(2), m.str(3), m.str(4), m.str(5)};
    }

    /* -- statePropertyDataType --
     * Deduces whether an entry in the stateProperties section describes loading
     * of state properties by direct value, file or function. The result is
     * returned as a StatePropertyDataType enumeration type.
     */

    static StatePropertyDataType statePropertyDataType(const std::string &propertyString, std::string &buffer)
    {
        static const std::regex reProperty(R"(.*=\s*(.+?)\s*$)");
        static const std::regex reValue(R"([\.0-9]+)");

        std::smatch m;

        if (std::regex_search(propertyString, m, reProperty))
        { // value or function
            buffer = m.str(1);

            if (std::regex_match(buffer, reValue))
            {
                return StatePropertyDataType::direct;
            }

            return StatePropertyDataType::function;
        }

        buffer = propertyString;

        return StatePropertyDataType::file;
    }

    /* -- propertyFunctionAndArguments --
     * Extracts the property function name and arguments from a given string and
     * stores them in two separate strings.
     */

    static bool propertyFunctionAndArguments(const std::string &totalString, std::string &functionName,
                                             std::string &argumentString)
    {
        static const std::regex reFuncArgs(R"(\s*(\w+)@?(.*))");
        std::smatch m;

        if (!std::regex_match(totalString, m, reFuncArgs))
            return false;

        functionName = m.str(1);
        argumentString = m.str(2);

        return true;
    }

    /* -- argumentsFromString --
     * Extracts the separate arguments from a string containing the arguments, finds
     * their corresponding double values and pushes them into the arguments vector. E.g.
     * "[gasTemperature, 1000]" first yields "gasTemperature" which is looked up in
     * the argumentMap to obtain its value, and then yields "1000" which is converted
     * into a double and pushed into the arguments vector.
     */

    static bool argumentsFromString(const std::string &argumentString, std::vector<double> &arguments,
                                    const std::map<std::string, double *> &argumentMap)
    {
        static const std::regex r(R"(\s*([\w\.]+)\s*(?:[,\]]|$))");

        for (auto it = std::sregex_iterator(argumentString.begin(), argumentString.end(), r);
             it != std::sregex_iterator(); ++it)
        {
            std::string current = it->str(1);

            if (isNumerical(current))
            {
                double value;

                if (!getValue(current, value))
                    return false;

                arguments.emplace_back(value);
            }
            else
            {
                if (argumentMap.count(current) == 0)
                    return false;

                arguments.emplace_back(*argumentMap.at(current));
            }
        }

        return true;
    }

    /* -- stateAndValue --
     * Extracts the state and the corresponding value from a line as provided in a
     * property file. They are both stored in their own string variables which are
     * passed by reference.
     */

    static bool stateAndValue(const std::string &propertyFileLine, std::string &stateString, std::string &valueString)
    {
        static const std::regex r(R"((.*?)\s*([\d\.e+-]+)\s*(?:\n|$))");
        std::smatch m;

        if (!std::regex_search(propertyFileLine, m, r))
            return false;

        stateString = m.str(1);
        valueString = m.str(2);

        return true;
    }

    /* -- statePropertyFile --
     * Parses a state property file into a vector of StateEntry, double pairs. This vector
     * is passed by reference.
     */

    static bool statePropertyFile(const std::string &fileName, std::vector<std::pair<StateEntry, double>> &entries)
    {
        const std::string inputPath = INPUT "/";

        std::ifstream in(inputPath + fileName);

        if (!in.is_open())
            return false;

        std::string line;

        while (std::getline(in, line))
        {
            line = removeComments(line);

            if (line.size() < 3)
                continue;

            std::string stateString, valueString;

            if (!Parse::stateAndValue(line, stateString, valueString))
                return false;

            double value;

            if (!Parse::getValue(valueString, value))
                return false;

            entries.emplace_back(Parse::propertyStateFromString(stateString), value);
        }

        return true;
    }

    /* -- rawCrossSectionFromStream --
     * Accepts a reference to an input file stream of an LXCat file. This stream should be
     * at a position just after reading a collision description from the LXCat file, since
     * this function searches for the line containing solely dashes, indicating that a
     * cross section follows. The raw cross section is then stored as a vector of pairs of
     * doubles, which the user passes by reference.
     */

    static void rawCrossSectionFromStream(std::vector<double> &rawEnergyData, std::vector<double> &rawCrossSection,
                                          std::istream &in)
    {
        std::string line;

        while (std::getline(in, line))
        {
            if (line.substr(0, 2) == "--")
                break;
        }

        double energy = 0., value = 0.;

        while (std::getline(in, line))
        {
            std::stringstream ss(line);
            if (!((ss >> energy) && (ss >> value)))
                break;

            rawEnergyData.emplace_back(energy);
            rawCrossSection.emplace_back(value);
        }
    }

    /* -- stringBufferFromFile --
     * Loads the complete content of a specified file into the given std::string.
     */

    static bool stringBufferFromFile(const std::string &fileName, std::string &buffer)
    {
        std::ifstream in(INPUT "/" + fileName);

        if (!in)
            return false;

        std::stringstream ss;
        ss << in.rdbuf();

        buffer = removeComments(ss.str());

        return true;
    }

    /* -- gasProperty --
     * Tries to parse the value of a given property for a given gas. The 'content' string
     * should store the content of the database file corresponding to the passed property
     * (i.e. mass or anharmonicFrequency).
     */

    static bool gasProperty(const std::string &gasName, double &property, const std::string &content)
    {
        const std::regex r(R"((?:^|\n))" + gasName + R"(\s+(.*)\s*)");
        std::smatch m;

        if (!std::regex_search(content, m, r))
            return false;

        std::stringstream ss(m[1]);
        return static_cast<bool>(ss >> property);
    }

    /* -- getValue --
     * Tries to parse a string into a double, returns a boolean based on its success
     * to do so.
     */

    static bool getValue(const std::string &valueString, double &value)
    {
        //            const std::regex r(R"(\s*(\d+\.?\d*)\s*\n*)");
        //            std::smatch m;
        //
        //            if (!std::regex_match(valueString, r)) return false;

        std::stringstream ss(valueString);

        return static_cast<bool>(ss >> value);
    }

    static bool isNumerical(const std::string &str)
    {
        static const std::regex reNum(R"(\s*\d*\.?\d+\s*)");

        return std::regex_match(str, reNum);
    }

    static Range getRange(const std::string &rangeString, bool &success)
    {
        static const std::regex r(
            R"(\s*((?:logspace\()|(?:linspace\())\s*(-?\d+\.?\d*)\s*,\s*(-?\d+\.?\d*)\s*,\s*(\d+\.?\d*))");
        std::smatch m;

        // No checking since this has already been performed once this function is called.
        if (!std::regex_search(rangeString, m, r))
        {
            success = false;
            return Range(0., 0., 0, false);
        }

        std::stringstream ss;

        const std::string function = m.str(1);

        double start, stop;
        uint32_t steps;

        ss << m[2];
        ss >> start;
        ss.clear();
        ss << m[3];
        ss >> stop;
        ss.clear();
        ss << m[4];
        ss >> steps;

        success = true;

        return Range(start, stop, steps, function[1] == 'o');
    }

  private:
    /* -- PARSING RANGES -- */

    /* -- firstValueInRange --
     * Tries to parse the first value in a string defining a range into a double,
     * returns a boolean based on its success to do so.
     */

    static bool firstValueInRange(const std::string &rangeString, double &value)
    {
        static const std::regex r(R"(\s*((?:logspace\(-?)|(?:linspace\())\s*(\d+\.?\d*)\s*)");
        std::smatch m;

        if (!std::regex_search(rangeString, m, r))
            return false;

        std::stringstream ss(m[2]);

        const std::string function = m.str(1);

        if (function[1] == 'o')
        {
            double power;
            bool success = static_cast<bool>(ss >> power);

            if (function.back() == '-')
                power = -power;

            value = std::pow(10., power);

            return success;
        }

        return static_cast<bool>(ss >> value);
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
template <>
inline bool Parse::setField<std::string>(const std::string &sectionContent, const std::string &fieldName,
                                         std::string &value)
{

    return getFieldValue(sectionContent, fieldName, value);
}

template <>
inline bool Parse::setField<bool>(const std::string &sectionContent, const std::string &fieldName, bool &value)
{
    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;

    std::stringstream s(valueBuffer);
    s >> std::boolalpha >> value;

    return true;
}

template <>
inline bool Parse::setField<Enumeration::EedfType>(const std::string &sectionContent, const std::string &fieldName,
                                                   Enumeration::EedfType &value)
{

    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;
    value = Enumeration::getEedfType(valueBuffer);
    return true;
}

template <>
inline bool Parse::setField<Enumeration::IonizationOperatorType>(const std::string &sectionContent,
                                                                 const std::string &fieldName,
                                                                 Enumeration::IonizationOperatorType &value)
{

    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;
    value = Enumeration::getIonizationOperatorType(valueBuffer);
    return true;
}

template <>
inline bool Parse::setField<Enumeration::GrowthModelType>(const std::string &sectionContent,
                                                          const std::string &fieldName,
                                                          Enumeration::GrowthModelType &value)
{

    std::string valueBuffer;

    if (!getFieldValue(sectionContent, fieldName, valueBuffer))
        return false;
    value = Enumeration::getGrowthModelType(valueBuffer);
    return true;
}

template <>
inline bool Parse::setField<std::vector<std::string>>(const std::string &sectionContent, const std::string &fieldName,
                                                      std::vector<std::string> &value)
{
    std::string fieldContent;

    if (!Parse::getSection(sectionContent, fieldName, fieldContent))
        return false;

    return Parse::getList(fieldContent, fieldName, value);
}

} // namespace loki

#endif // LOKI_CPP_PARSE_H
