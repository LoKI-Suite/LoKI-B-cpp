//
// Created by daan on 2-5-19.
//

#include <iostream>

#include <Enumeration.h>
#include <Constant.h>
#include <Setup.h>

using std::cout;
using std::endl;

using namespace loki::Enumeration;

int main (int argc, char ** argv)
{
    // TODO: Allow the user to specify the setup file in the program arguments.;

    loki::Setup setup;

    if (!setup.parseFile("default_lokib_setup.in")) {
        return 1;
    }

    if (setup.electronKinetics.isOn) {
        if (setup.electronKinetics.eedfType == EedfType::boltzmann) {
            std::cout << "solver routine" << std::endl;
        }
    }

    // generate output

    return 0;
}