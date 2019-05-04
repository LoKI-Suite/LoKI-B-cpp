//
// Created by daan on 2-5-19.
//

#include <iostream>
#include "Setup.h"

namespace loki {
    using namespace Enumeration;

    Setup::Setup(const std::string& fileName) {
        std::cout << "Parsing setup file: " << fileName << std::endl;
    }
}