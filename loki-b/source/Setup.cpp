//
// Created by daan on 2-5-19.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "Setup.h"

namespace loki {
    Setup::Setup(const std::string& fileName) {
        std::cout << "Parsing setup file: " << inputPath + '/' + fileName << std::endl;

        std::ifstream file(inputPath + '/' + fileName);

        if (!file.is_open()) {
            std::cerr << "could not open input file" << std::endl;
            return;
        }

        std::stringstream stringBuffer;
        stringBuffer << file.rdbuf();

        std::string fileContent = stringBuffer.str();

        // This regular expression matches a specific section in the input file,
        // e.g. the "electronKinetics" section.
        std::regex r(R"(electronKinetics:\s*\n([^]*?)((\n+\w)|(\n+$)))");
        std::smatch m;

        if(std::regex_search(fileContent, m, r)) {

            std::cout << m[1] << std::endl;

        } else {

            std::cout << "could not match regexp" << std::endl;

        }
    }
}