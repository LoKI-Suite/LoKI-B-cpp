//
// Created by daan on 2-5-19.
//

#include <iostream>

#include <Enumeration.h>
#include <Constant.h>
#include <Setup.h>
#include <Simulation.h>
#include <chrono>
#include <Log.h>

using std::cout;
using std::endl;

using namespace loki::Enumeration;

int main (int argc, char ** argv)
{
    // TODO: Allow the user to specify the setup file in the program arguments.

    auto begin = std::chrono::high_resolution_clock::now();

    try {
        loki::Setup setup;

        if (!setup.parseFile("default_lokib_setup.in")) {
            return 1;
        }

        loki::Simulation simulation(setup);

        if (setup.electronKinetics.isOn) {
            if (setup.electronKinetics.eedfType == EedfType::boltzmann) {
                std::cout << "solver routine" << std::endl;
            }
        }
    } catch (const std::exception &e) {
        std::cout << e.what() << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << "mus" << std::endl;

    // generate output

    return 0;
}