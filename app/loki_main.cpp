/** \file
 *
 *  Implementation of the LoKI-B simulation program.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2022 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   2 May 2019 (first C++ version)
 */

#include "LoKI-B/LegacyToJSON.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Simulation.h"
#include "LoKI-B/Gnuplot.h"
#include "LoKI-B/Output.h"
#include <chrono>
#include <exception>
#include <iostream>

//#define LOKIB_ENABLE_FPU_EXCEPTIONS

#ifdef LOKIB_ENABLE_FPU_EXCEPTIONS
#include <cfenv>
#endif

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::WorkingConditions &wc,
                   const loki::Power &power, const loki::EedfCollisionDataMixture &coll_data,
                   const loki::SwarmParameters &swarmParameters,
                   const loki::Vector *firstAnisotropy)
{
    using namespace loki;
    writeGnuplot(std::cout, "Eedf", "Energy (eV)", "Eedf (Au)", grid.getCells(), eedf);
    // writeGnuplot("Cross Section", "Energy (eV)", "Cross Section (m^{-2})", grid.getNodes(),
    // *gases[0]->collisions[(uint8_t)loki::Enumeration::CollisionType::excitation][0]->crossSection);

    // loki::Vector *nodeCrossSection =
    // gases[0]->collisions[(uint8_t)loki::Enumeration::CollisionType::excitation][0]->crossSection;

    // loki::Vector cellCrossSection(grid.nCells());

    // for (loki::Grid::Index i = 0; i < grid.nCells(); ++i)
    //     cellCrossSection[i] = 0.5 * ((*nodeCrossSection)[i] + (*nodeCrossSection)[i + 1]);

    // writeGnuplot("Cell Cross Section", "Energy (eV)", "Cross Section (m^{-2})", grid.getCells(), cellCrossSection);

    // for (loki::Grid::Index i = 0; i < grid.nCells(); ++i)
    // {
    //     if (cellCrossSection[i] != 0.)
    //         std::cerr << grid.getCell(i) << "\t" << cellCrossSection[i] << '\n';
    // }

    // if (firstAnisotropy)   writeGnuplot("First Anisotropy", "Energy (eV)", "First Anisotropy (Au)", grid.getCells(), firstAnisotropy);
}

void handleExistingOutputPath(std::string &folder)
{
    std::cerr << "Please enter a new destination for the output files (keep empty for unaltered).\nOutput/";
    std::getline(std::cin, folder);
}

const char* usage_str = "Usage: loki [--convert] <inputfile>";

int main(int argc, char **argv)
try
{
#ifdef LOKIB_ENABLE_FPU_EXCEPTIONS
    /* see https://en.cppreference.com/w/cpp/numeric/fenv for making this portable. There are
     * also other possibilities, like translating the FPU exception into a C++ exception that
     * can then be handled gracefully by the program... so it seems.
     * NOTE that division by zero (of a non-zero value) and overflows are not necessarily
     *      problematic, since IEE754 describes clear semantics for these cases.
     */
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    /* The target for JSON output. This will be initialized if JSON output
     * is asked for (see the cnf.at("output") bit below). In that case this
     * object will be written to console at the end of the simulation (not
     * yet to file). Note that support is incomplete, see Output.cpp.
     */
    std::unique_ptr<loki::json_type> json_output_data;

    auto begin = std::chrono::high_resolution_clock::now();
    bool convert_input = false;
    std::filesystem::path fileName;
    if (argc==2)
    {
            convert_input = false;
            fileName = argv[1];
    }
    else if (argc==3)
    {
            if (argv[1]!=std::string{"--convert"})
            {
                throw std::runtime_error(usage_str);
            }
            convert_input = true;
            fileName = argv[2];
    }
    else
    {
        throw std::runtime_error(usage_str);
    }

    std::vector<std::unique_ptr<loki::Output>> output;
    /* Arguments: [--convert] filename
     * - If filename has extension '.in', conversion to JSON is always done and
     *   option --convert is ignored when specified.
     * - If filename has extensio .json, --convert can be used to patch a 'legacy'
     *   JSON file (a 'literal' translation of the '.in' file to JSON) to bring it
     *   in the required form. When the JSON file already is already in the 'new'
     *   JSON format, do NOT pass --convert.
     */
    const bool input_is_json = fileName.has_extension() && fileName.extension() == ".json";
    const loki::json_type cnf = input_is_json
        ? (convert_input ? loki::legacyToJSON(loki::read_json_from_file(fileName)) : loki::read_json_from_file(fileName))
        : loki::legacyToJSON(fileName);
    if (convert_input)
    {
        if (!input_is_json)
        {
            std::cout << "Warning: option --convert is ignored: conversion is "
                         "implicit for legacy ('.in') input files." << std::endl;
        }
        else
        {
            std::cout << "Input file << '" << fileName << "' was converted. "
                         "Result:\n" << cnf.dump(2) << std::endl;
        }
    }
    loki::Simulation simulation(fileName, cnf);
    if (cnf.at("output").at("isOn"))
    {
        // Write text files by default: in the absence of "writeText" or
        // when such key exists and its value is true.
        if (cnf.at("output").contains("writeText")==false || cnf.at("output").at("writeText"))
        {
            output.emplace_back(
                new loki::FileOutput(cnf, &simulation.workingConditions(),
                        &handleExistingOutputPath));
        }
        // Write JSON (to the console, at the end) only when asked for: when key 'writeJSON' exists
        // AND the value is true.
        if (cnf.at("output").contains("writeJSON")==true && cnf.at("output").at("writeJSON"))
        {
            json_output_data.reset(new loki::json_type);
            output.emplace_back(
                new loki::JsonOutput(*json_output_data, cnf, &simulation.workingConditions()));
        }
    }
    // register all output producers with the simulation's obtainedResults event
    for (const auto& out : output)
    {
        simulation.obtainedResults().addListener(&loki::Output::saveCycle, out.get());
    }

    simulation.obtainedResults().addListener(handleResults);

    simulation.run();
    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "Simulation finished, elapsed time = "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
              << "mus" << std::endl;

    // if data wer harvested (also) in JSON form, print the JSON output object to screen.
    if (json_output_data)
    {
        std::cout << "Output data (JSON format):" << std::endl;
        std::cout << json_output_data->dump(2) << std::endl;
    }

    return 0;
}
catch (const std::exception &e)
{
    std::cerr << e.what() << std::endl;
    return 1;
}
