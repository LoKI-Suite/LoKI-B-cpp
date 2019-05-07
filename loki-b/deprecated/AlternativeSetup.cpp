//
// Created by daan on 6-5-19.
//

#include "AlternativeSetup.h"

#include <regex>
#include <iostream>
#include <fstream>

#include "Parse.h"

namespace loki {
    template <> std::string &SetupStruct::add<std::string>(const std::string &name, const std::string &value) {
        vString.emplace_back(value);
        key.emplace(name, vString.size()-1);

        return vString.back();
    }
    template <> SetupStruct &SetupStruct::add<SetupStruct>(const std::string &name, const SetupStruct &value) {
        vSetupStruct.emplace_back(value);
        key.emplace(name, vSetupStruct.size()-1);

        return vSetupStruct.back();
    }

    template <>
    std::string &SetupStruct::get<std::string>(const std::string &name) {
        const auto &id = key[name];
        return vString[id];
    }

    template <>
    SetupStruct &SetupStruct::get<SetupStruct>(const std::string &name) {
        const auto &id = key[name];
        return vSetupStruct[id];
    }

    SetupStruct &SetupStruct::operator[](const std::string &name) {
        return get<SetupStruct>(name);
    }

    std::string &SetupStruct::operator()(const std::string &name) {
        return get<std::string>(name);
    }

    bool SetupStruct::parseFile(const std::string &fileName) {
        std::ifstream file(inputPath + '/' + fileName);

        if (!file.is_open()) {
            std::cerr << "Could not find/open specified file" << std::endl;
            return false;
        }

        // Load the contents of the file into a string buffer.
        std::stringstream stringBuffer;
        stringBuffer << file.rdbuf();

        // Store the file contents in a string and remove any comments.
        std::string fileContent = Parse::removeComments(stringBuffer.str());

        return parseSection(fileContent);
    }

    bool SetupStruct::parseSection(const std::string &content) {
        const std::string sLevelString = R"((?:[ ]{)" + std::to_string(level*2) + R"(}))";
        const std::string eLevelString = R"((?:[ ]{0,)" + std::to_string(level*2) + R"(}))";

        const std::regex r(sLevelString + R"((\w+):\s*\n*([^]*?)(?:(?:\n*\s*$)|(?:\n+)" + eLevelString + R"([a-zA-Z])))");
        std::smatch ma;

        std::string spaces;

        for (uint8_t i = 0; i < level+1; ++i) spaces += "  ";

        uint32_t index = 0;

        if (std::regex_search(content.begin(), content.end(), ma, r)) {
            do {
                std::string sectionContent = spaces;
                sectionContent += ma[2];

                if (sectionContent.find('\n') != std::string::npos) {
                    // Content is multiline.
//                    std::cout << spaces << ma[1] << ", " << level+1 << std::endl;
                    if (!add(ma[1], SetupStruct(level+1)).parseSection(sectionContent))
                        return false;
                } else {
                    // Content is a value (base case).
                    if (sectionContent.empty()) return false;

//                    std::cout << spaces << ma[1] << "\tbase" << std::endl;

                    add(ma[1], sectionContent);
                }

                index += ma.position(2) + ma[2].length();

            } while (std::regex_search(content.begin() + index, content.end(), ma, r));
        } else {
            if (content.empty()) return false;

            const std::regex re(sLevelString + R"(-\s*([^]*?)(?:(?:\n*\s*$)|(?:\n+)" + eLevelString + R"([-a-zA-Z])))");

            while (std::regex_search(content.begin() + index, content.end(), ma, re)) {

                const std::string value = ma[1];

                add(std::to_string(index), value);

//                std::cout << spaces << value << std::endl;

                index += ma.position(1) + ma[1].length();
            }
        }

        return true;
    }
}