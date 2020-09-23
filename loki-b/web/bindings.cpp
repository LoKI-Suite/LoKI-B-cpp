#include <iostream>

#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/Simulation.h"
#include <chrono>
#include <exception>

#include <emscripten/bind.h>

using loki::Log;
using loki::Message;

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::WorkingConditions &wc,
                   const loki::Power &power, const std::vector<loki::EedfGas *> &gases,
                   const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients,
                   const loki::Vector &firstAnisotropy);

void handleExistingOutputPath(std::string &folder);

void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel, const loki::Vector &x,
          const loki::Vector &y);

int run(std::string file_name)
{
    try
    {
        auto begin = std::chrono::high_resolution_clock::now();

        std::unique_ptr<loki::Simulation> simulation;
        if (file_name.size() >= 5 && file_name.substr(file_name.size() - 5) == ".json")
        {
            const loki::json_type cnf = loki::read_json_from_file(file_name);
            simulation.reset(new loki::Simulation(cnf));
        }
        else
        {
            return 1;
        }

        simulation->obtainedResults.addListener(handleResults);
        simulation->outputPathExists.addListener(handleExistingOutputPath);

        simulation->run();
        auto end = std::chrono::high_resolution_clock::now();
        Log<Message>::Notify("Simulation finished, elapsed time = ",
                             std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count(), "mus");

        // generate output

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}

EMSCRIPTEN_BINDINGS(loki)
{
    emscripten::function("run", &run);
}

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::WorkingConditions &wc,
                   const loki::Power &power, const std::vector<loki::EedfGas *> &gases,
                   const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients, const loki::Vector &firstAnisotropy)
{
    plot("Eedf", "Energy (eV)", "Eedf (Au)", grid.getCells(), eedf);
}

void handleExistingOutputPath(std::string &folder)
{
    std::cerr << "Please enter a new destination for the output files (keep empty for unaltered).\nOutput/";
    std::getline(std::cin, folder);
}

void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel, const loki::Vector &x,
          const loki::Vector &y)
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
