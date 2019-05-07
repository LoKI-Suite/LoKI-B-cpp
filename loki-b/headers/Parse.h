//
// Created by daan on 6-5-19.
//

#ifndef LOKI_CPP_PARSE_H
#define LOKI_CPP_PARSE_H

#include <string>
#include <regex>

#include "Enumeration.h"

namespace loki {
    struct Parse {

        // TODO: add commenting for the setField function.

        template <typename T>
        static bool setField(const std::string &sectionContent,
                             const std::string &fieldName, T &value) {

            std::string valueBuffer;

            if (!getFieldValue(sectionContent, fieldName, valueBuffer))
                return false;

            std::stringstream s(valueBuffer);
            s >> value;

            return true;
        }

        // TODO: add commenting for the getFieldValue function.

        static bool getFieldValue(const std::string &sectionContent,
                const std::string &fieldName, std::string &valueBuffer) {

            std::regex r(fieldName + R"(:\s*([^\s\n]*)\s*\n*)");
            std::smatch m;

            if (!std::regex_search(sectionContent, m, r))
                return false;

            valueBuffer = m[1];

            return true;
        }

        // TODO: add commenting for the getList function.

        static bool getList(const std::string &sectionContent,
                const std::string &fieldName, std::vector<std::string> &container) {

            std::regex r(R"(-\s*(.*)\s*\n*)");

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
         *
         * NOTE: This function does not deal with ambiguous section names. E.g. the "isOn" field
         * occurs multiple times in the file. The user is advised to only retrieve sections that
         * are one level above the level of "fileContent". E.g. retrieve "isOn" when "fileContent"
         * contains the contents of the "electronKinetics" section.
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
            const std::regex reSection(sectionTitle + R"(:\s*\n*([^]*?)(?:(?:\n+)" + levelString + R"(\w)|(?:\n*\s*$)))");

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

    };

    // TODO: Add commenting for the specialized setField functions.

    // Use "inline" with fully specialized templates to refrain from
    // breaking the "One Definition Rule".
    template <> inline
    bool Parse::setField<std::string>(const std::string &sectionContent,
                               const std::string &fieldName, std::string &value) {

        return getFieldValue(sectionContent, fieldName, value);
    }

    template <> inline
    bool Parse::setField<bool>(const std::string &sectionContent,
                               const std::string &fieldName, bool &value) {
        std::string valueBuffer;

        if (!getFieldValue(sectionContent, fieldName, valueBuffer))
            return false;

        std::stringstream s(valueBuffer);
        s >> std::boolalpha >> value;

        return true;
    }

    template <> inline
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

    template <> inline
    bool Parse::setField<Enumeration::IonizationOperatorType>(const std::string &sectionContent,
            const std::string &fieldName, Enumeration::IonizationOperatorType &value) {

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

    template <> inline
    bool Parse::setField<Enumeration::GrowthModelType>(const std::string &sectionContent,
            const std::string &fieldName, Enumeration::GrowthModelType &value) {

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

    template <> inline
    bool Parse::setField<std::vector<std::string>>(const std::string &sectionContent,
            const std::string &fieldName, std::vector<std::string> &value) {

        std::string fieldContent;

        if (!Parse::getSection(sectionContent, fieldName, fieldContent))
            return false;

        return Parse::getList(fieldContent, fieldName, value);

    }
}

#endif //LOKI_CPP_PARSE_H
