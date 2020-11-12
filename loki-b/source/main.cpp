//
// Created by daan on 2-5-19.
//

#include <iostream>

#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/Simulation.h"
#include "LoKI-B/Gnuplot.h"
#include <chrono>
#include <exception>

//#define LOKIB_ENABLE_FPU_EXCEPTIONS

#ifdef LOKIB_ENABLE_FPU_EXCEPTIONS
#include <cfenv>
#endif

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
                   const loki::Power &power, const std::vector<loki::EedfGas *> &gases,
                   const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients, const loki::Vector &firstAnisotropy)
{
    using namespace loki;
    writeGnuplot(std::cout, "Eedf", "Energy (eV)", "Eedf (Au)", grid.getCells(), eedf);
    // writeGnuplot("Cross Section", "Energy (eV)", "Cross Section (m^{-2})", grid.getNodes(),
    // *gases[0]->collisions[(uint8_t)loki::Enumeration::CollisionType::excitation][0]->crossSection);

    // loki::Vector *nodeCrossSection =
    // gases[0]->collisions[(uint8_t)loki::Enumeration::CollisionType::excitation][0]->crossSection;

    // loki::Vector cellCrossSection(grid.nCells());

    // for (uint32_t i = 0; i < grid.nCells(); ++i)
    //     cellCrossSection[i] = 0.5 * ((*nodeCrossSection)[i] + (*nodeCrossSection)[i + 1]);

    // writeGnuplot("Cell Cross Section", "Energy (eV)", "Cross Section (m^{-2})", grid.getCells(), cellCrossSection);

    // for (uint32_t i = 0; i < grid.nCells(); ++i)
    // {
    //     if (cellCrossSection[i] != 0.)
    //         std::cerr << grid.getCell(i) << "\t" << cellCrossSection[i] << '\n';
    // }

    //    writeGnuplot("First Anisotropy", "Energy (eV)", "First Anisotropy (Au)", grid.getCells(), firstAnisotropy);
}

void handleExistingOutputPath(std::string &folder)
{
    std::cerr << "Please enter a new destination for the output files (keep empty for unaltered).\nOutput/";
    std::getline(std::cin, folder);
    //    std::cin >> folder;
}


int main(int argc, char **argv)
{
#ifdef LOKIB_ENABLE_FPU_EXCEPTIONS
    // see https://en.cppreference.com/w/cpp/numeric/fenv for making this portable. There are
    // also other possibilities, like translating the FPU exception into a C++ exception that
    // can then be handled gracefully by the program... so it seems.
    // NOTE that division by zero (of a non-zero value) and overflows are not necessarily
    //      problematic, since IEE754 describes clear semantics for thee cases.
    std::feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    try
    {
        /* When WRITE_OUTPUT_TO_JSON_OBJECT is defined, a JSONOutput object will
         * be set up instead of af FileOutput object. The variable data_out will
         * act as its output root object. Note that support is incomplete, see
         * Output.cpp.
         */
//#define WRITE_OUTPUT_TO_JSON_OBJECT
#ifdef WRITE_OUTPUT_TO_JSON_OBJECT
        loki::json_type data_out;
#endif
        auto begin = std::chrono::high_resolution_clock::now();
        if (argc != 2)
        {
            throw std::runtime_error("Expected the input file as the single argument.");
        }

        std::unique_ptr<loki::Simulation> simulation;
        std::string fileName(argv[1]);
        if (fileName.size() >= 5 && fileName.substr(fileName.size() - 5) == ".json")
        {
            const loki::json_type cnf = loki::read_json_from_file(fileName);
            simulation.reset(new loki::Simulation(cnf));
            if (cnf.at("output").at("isOn"))
            {
                loki::Output* output =
#ifdef WRITE_OUTPUT_TO_JSON_OBJECT
                    new loki::JsonOutput(data_out, cnf, &simulation->m_workingConditions,
                            &simulation->m_jobManager);
#else
                    new loki::FileOutput(cnf, &simulation->m_workingConditions,
                            &simulation->m_jobManager, &handleExistingOutputPath);
#endif
                simulation->configureOutput(output);
            }
        }
        else
        {
            const loki::Setup setup(argv[1]);
            simulation.reset(new loki::Simulation(setup));
            if (setup.output.isOn)
            {
                loki::Output* output = new loki::FileOutput(setup, &simulation->m_workingConditions,
                        &simulation->m_jobManager, &handleExistingOutputPath);
                simulation->configureOutput(output);
            }
        }

        simulation->m_obtainedResults.addListener(handleResults);

        simulation->run();
        auto end = std::chrono::high_resolution_clock::now();
        std::cerr << "Simulation finished, elapsed time = "
                  << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                  << "mus" << std::endl;

#ifdef WRITE_OUTPUT_TO_JSON_OBJECT
        std::cout << "Output data:" << std::endl;
        std::cout << data_out.dump(2) << std::endl;
#endif

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}

