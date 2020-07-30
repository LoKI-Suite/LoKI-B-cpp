//
// Created by daan on 2-5-19.
//

#include <iostream>

#include <Setup.h>
#include <Simulation.h>
#include <LinearAlgebra.h>
#include <chrono>
#include <exception>

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

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::WorkingConditions &wc,
                   const loki::Power &power, const std::vector<loki::EedfGas*>& gases,
                   const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients,
                   const loki::Vector &firstAnisotropy);

void handleExistingOutputPath(std::string &folder);

void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
          const loki::Vector &x, const loki::Vector &y);

int main(int argc, char **argv)
{
    try
    {
        auto begin = std::chrono::high_resolution_clock::now();
        if (argc != 2)
        {
            throw std::runtime_error("Expected the input file as the single argument.");
        }

        std::unique_ptr<loki::Simulation> simulation;
        std::string fileName(argv[1]);
        if (fileName.size()>=5 && fileName.substr(fileName.size()-5)==".json")
        {
            const loki::json_type cnf = loki::read_json_from_file(fileName);
            simulation.reset(new loki::Simulation(cnf));
        }
        else
        {
            const loki::Setup setup(argv[1]);
            simulation.reset(new loki::Simulation(setup));
        }

        simulation->obtainedResults.addListener(handleResults);
        simulation->outputPathExists.addListener(handleExistingOutputPath);

        simulation->run();
        auto end = std::chrono::high_resolution_clock::now();
        std::cerr << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;

        // generate output

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::WorkingConditions &wc,
                   const loki::Power &power, const std::vector<loki::EedfGas*>& gases,
                   const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients,
                   const loki::Vector &firstAnisotropy)
{
    plot("Eedf", "Energy (eV)", "Eedf (Au)", grid.getCells(), eedf);
    // plot("Cross Section", "Energy (eV)", "Cross Section (m^{-2})", grid.getNodes(), *gases[0]->collisions[(uint8_t)loki::Enumeration::CollisionType::excitation][0]->crossSection);

    // loki::Vector *nodeCrossSection = gases[0]->collisions[(uint8_t)loki::Enumeration::CollisionType::excitation][0]->crossSection;

    // loki::Vector cellCrossSection(grid.cellNumber);

    // for (uint32_t i = 0; i < grid.cellNumber; ++i)
    //     cellCrossSection[i] = 0.5 * ((*nodeCrossSection)[i] + (*nodeCrossSection)[i + 1]);

    // plot("Cell Cross Section", "Energy (eV)", "Cross Section (m^{-2})", grid.getCells(), cellCrossSection);

    // for (uint32_t i = 0; i < grid.cellNumber; ++i)
    // {
    //     if (cellCrossSection[i] != 0.)
    //         std::cerr << grid.getCell(i) << "\t" << cellCrossSection[i] << '\n';
    // }

    //    this->plot("First Anisotropy", "Energy (eV)", "First Anisotropy (Au)", grid.getCells(), firstAnisotropy);
}

void handleExistingOutputPath(std::string &folder)
{
    std::cerr << "Please enter a new destination for the output files (keep empty for unaltered).\nOutput/";
    std::getline(std::cin, folder);
    //    std::cin >> folder;
}

void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
          const loki::Vector &x, const loki::Vector &y)
{
    std::cout << "unset key" << std::endl;
    std::cout << "set xlabel \"" << xlabel << "\"" << std::endl;
    std::cout << "set ylabel \"" << ylabel << "\"" << std::endl;
    std::cout << "set title \"" << title << "\"" << std::endl;
    std::cout << "set xrange [" << x[0] << ":" << x[x.size() - 1] << "]" << std::endl;
    std::cout << "set logscale y" << std::endl;
    // std::cout << "set format y '%g'" << std::endl;
    std::cout << "plot '-' w l" << std::endl;
    for (uint32_t i = 0; i < x.size(); ++i)
    {
        std::cout << x[i] << "\t" << y[i] << '\n';
    }
    std::cout << "e" << std::endl;
}
