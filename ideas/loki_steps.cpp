/** \file
 *
 *  Implementation of the LoKI-B simulation program.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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

#include <iostream>

#include "LoKI-B/LegacyToJSON.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Simulation.h"
#include "LoKI-B/Gnuplot.h"
#include "LoKI-B/Output.h"
#include <chrono>
#include <exception>
#include <fstream>

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
    //    std::cin >> folder;
}


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

    auto begin = std::chrono::high_resolution_clock::now();
    if (argc != 2)
    {
        throw std::runtime_error("Usage: loki <inputfile>");
    }
    std::filesystem::path fileName(argv[1]);
    const bool input_is_json = fileName.has_extension() && fileName.extension() == ".json";
    const loki::json_type cnf = input_is_json
        ? loki::read_json_from_file(fileName)
        : loki::legacyToJSON(fileName);

    loki::WorkingConditions workingConditions(cnf.at("workingConditions"));
    std::unique_ptr<loki::ElectronKinetics> electron_kinetics;
    const std::string eedfType = cnf.at("electronKinetics").at("eedfType");
    if (eedfType=="boltzmann")
    {
	electron_kinetics.reset(new loki::ElectronKineticsBoltzmann(fileName, cnf.at("electronKinetics"), &workingConditions));
    }
    else if (eedfType=="prescribed")
    {
	electron_kinetics.reset(new loki::ElectronKineticsPrescribed(fileName, cnf.at("electronKinetics"), &workingConditions));
    }
    else
    {
        throw std::runtime_error("Bad value of 'eedfType': expected 'boltzmann' or 'prescribed'.");
    }

    std::unique_ptr<loki::Output> output;
    std::unique_ptr<loki::json_type> data_out;
    if (cnf.at("output").at("isOn"))
    {
        data_out.reset(new loki::json_type);
        output.reset(new loki::JsonOutput(*data_out, cnf, &workingConditions));
        electron_kinetics->obtainedNewEedf.addListener(&loki::Output::saveCycle, output.get());
    }
    electron_kinetics->obtainedNewEedf.addListener(handleResults);

    /* we bypass the job manager. Instead we select some case parameter values
     * (in this case only for E/N), install those values and call solve()
     * ourselves.
     */
    if (eedfType=="boltzmann")
    {
        for (double E_N : {0.1, 1., 10.})
        {
            std::cout << "Running loki for E/N = " << E_N << " Td" << std::endl;
            workingConditions.updateReducedField(E_N);
            std::stringstream id;
            id << "ReducedField_" << E_N;
            workingConditions.setCurrentJobFolder(id.str());
            electron_kinetics->solve();
        }
    }
    else if (eedfType=="prescribed")
    {
        for (double Te : {0.1, 1., 10.})
        {
            std::cout << "Running loki for Te = " << Te << " eV" << std::endl;
            workingConditions.updateElectronTemperature(Te);
            std::stringstream id;
            id << "Te_" << Te;
            workingConditions.setCurrentJobFolder(id.str());
            electron_kinetics->solve();
        }
    }
    else
    {
        throw std::runtime_error("Bad value of 'eedfType': expected 'boltzmann' or 'prescribed'.");
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "Simulation finished, elapsed time = "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
              << "mus" << std::endl;

    if (data_out.get())
    {
        std::string fname = "results.json";
        std::cout << "Output data (JSON format) will be written to file '"
            << fname << "'." << std::endl;
        std::ofstream ofs(fname);
        ofs << data_out->dump(2) << std::endl;
    }

    return 0;
}
catch (const std::exception &e)
{
    std::cerr << e.what() << std::endl;
    return 1;
}
