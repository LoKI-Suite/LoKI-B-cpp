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

// TODO: Cleanup
//  [DONE] 1. Check the equal sharing and one takes all ionization routines
//  [DONE] 2. Check attachment (with oxygen)
//  3. Switch from for loops to vector expressions
//  [DONE] 4. Fix bad_alloc when collisions have thresholds above umax (check before adding the collision).

// TODO: Steps to finishing the project:
//  1. [DONE] Write functions to compute swarm parameters and rate coefficients
//  2. Write a class to write the output to file
//  3. Implement the smart grid
//  4. Implement a simple version of the jobs
//  5. Separate backend and front end

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

        simulation.run();

    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count() << "mus" << std::endl;

    // generate output

    return 0;
}
