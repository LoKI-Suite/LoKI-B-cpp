//
// Created by daan on 2-5-19.
//

#include <iostream>

#include <Setup.h>
#include <Simulation.h>
#include <LinearAlgebra.h>
#include <chrono>
#include <Log.h>

// TODO: Cleanup
//  [DONE] 1. Check the equal sharing and one takes all ionization routines
//  [DONE] 2. Check attachment (with oxygen)
//  3. Switch from for loops to vector expressions
//  [DONE] 4. Fix bad_alloc when collisions have thresholds above umax (check before adding the collision).

// DONE: Steps to finishing the project:
//  1. [DONE] Write functions to compute swarm parameters and rate coefficients
//  2. [DONE] Write a class to write the output to file
//  3. [DONE] Implement the smart grid
//  4. [DONE] Implement jobs
//  5. [DONE] Separate backend and front end

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::Power &power,
                   const std::vector<loki::EedfGas *> &gasses, const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients,
                   const loki::Vector &firstAnisotropy);

void handleExistingOutputPath(std::string &folder);

void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
          const loki::Vector &x, const loki::Vector &y);

int main(int argc, char **argv) {
    if (argc != 2) {
        loki::Log<loki::Message>::Warning("Expected the input file as the single argument.");
        return 1;
    }

    std::string fileName = argv[1];

    auto begin = std::chrono::high_resolution_clock::now();

    try {
        loki::Setup setup;

        if (!setup.parseFile(fileName)) {
            return 1;
        }

        loki::Simulation simulation(setup);

        simulation.obtainedResults.addListener(handleResults);
        simulation.outputPathExists.addListener(handleExistingOutputPath);

        simulation.run();

    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;

    // generate output

    return 0;
}

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::Power &power,
                   const std::vector<loki::EedfGas *> &gasses, const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients,
                   const loki::Vector &firstAnisotropy) {
    plot("Eedf", "Energy (eV)", "Eedf (Au)", grid.getCells(), eedf);
//    this->plot("First Anisotropy", "Energy (eV)", "First Anisotropy (Au)", grid.getCells(), firstAnisotropy);
}

void handleExistingOutputPath(std::string &folder) {
    std::cerr << "Please enter a new destination for the output files (keep empty for unaltered).\nOutput/";
    std::getline(std::cin, folder);
//    std::cin >> folder;
}

void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
          const loki::Vector &x, const loki::Vector &y) {
    std::cout << "unset key" << std::endl;
    std::cout << "set xlabel \"" << xlabel << "\"" << std::endl;
    std::cout << "set ylabel \"" << ylabel << "\"" << std::endl;
    std::cout << "set title \"" << title << "\"" << std::endl;
    std::cout << "set xrange [" << x[0] << ":" << x[x.size() - 1] << "]" << std::endl;
    std::cout << "set logscale y" << std::endl;
    std::cout << "plot '-' w l" << std::endl;
    for (uint32_t i = 0; i < x.size(); ++i) {
        std::cout << x[i] << "\t" << y[i] << '\n';
    }
    std::cout << "e" << std::endl;
}
