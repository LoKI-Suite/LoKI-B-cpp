//
// Created by daan on 2-5-19.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "Setup.h"

namespace loki {
    // TODO: move the functionality that is now in the constructor into a Setup::parse
    //  function such that it can return a boolean to specify whether the operation was
    //  successful.

    // TODO: add "parse" functions for the substructures of the setup structure, such
    //  that the substructure can be filled by passing a section string.

    Setup::Setup(const std::string& fileName) {
        std::ifstream file(inputPath + '/' + fileName);

        if (!file.is_open()) {
            std::cerr << "could not open input file" << std::endl;
            return;
        }

        // Load the contents of the file into a string buffer.
        std::stringstream stringBuffer;
        stringBuffer << file.rdbuf();

        // Store the file contents in a string and remove any comments.
        std::string fileContent = Setup::removeComments(stringBuffer.str());

        // Object to store the results of the Setup:getSection calls.
        std::string sectionBuffer;

        if (Setup::getSection(fileContent, "workingConditions", sectionBuffer)) {
            //parse workingConditions section
            std::cout << sectionBuffer << std::endl;
        } else {
            std::cerr << "The input file does not contain a section to specify the working "
                         "conditions. Please add this section and try again." << std::endl;
            return;
        }

        if (Setup::getSection(fileContent, "electronKinetics", sectionBuffer)) {
            //parse electronKinetics section
            std::cout << sectionBuffer << std::endl;
        } else {
            this->electronKinetics.isEnabled = false;
        }

        if (Setup::getSection(fileContent, "output", sectionBuffer)) {
            //parse output section
            std::cout << sectionBuffer << std::endl;
        } else {
            std::cerr << "The input file does not contain a section to specify the output. "
                         "Please add this section and try again." << std::endl;
            return;
        }
    }

    bool Setup::getSection(const std::string &fileContent, const std::string &sectionTitle,
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
        const std::regex reSection(sectionTitle + R"(:\s*\n*([^]*?)\n+(?:(?:)" + levelString + R"(\w)|$))");

        if (!std::regex_search(fileContent, m, reSection))
            return false;

        sectionBuffer = levelString + "  ";
        sectionBuffer += m[1];

        return true;
    }

    std::string Setup::removeComments(const std::string &content) {
        const std::regex reClean(R"(%[^]*?\n)");
        return std::regex_replace(content, reClean, "\n");
    }
}