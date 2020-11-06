#include <iostream>

#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Setup.h"
#include "LoKI-B/Simulation.h"
#include "emscripten/em_asm.h"
#include <chrono>
#include <exception>
#include <sstream>

#include <emscripten.h>
#include <emscripten/bind.h>

/// \todo Explain why it is not possible to just use smth. like loki::json_type in the code.
using loki::Log;
using loki::Message;
using loki::json_type;

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::WorkingConditions &wc,
                   const loki::Power &power, const std::vector<loki::EedfGas *> &gases,
                   const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients,
                   const loki::Vector &firstAnisotropy)
{
    EM_ASM({ plot($0, $1, $2, $3); }, grid.getCells().data(), grid.getCells().size(), eedf.data(), eedf.size());
}

void handleExistingOutputPath(std::string &folder)
{
    std::cerr << "Please enter a new destination for the output files (keep empty for unaltered).\nOutput/";
    std::getline(std::cin, folder);
}

int run(std::string file_contents)
{
    try
    {
        auto begin = std::chrono::high_resolution_clock::now();

        std::unique_ptr<loki::Simulation> simulation;

        std::stringstream ss(file_contents);
        const loki::json_type cnf = loki::read_json_from_stream(ss);

        /* when non-null, a JSonOutput object will be created to produce output,
         * rather than a FileOutput object. For the web version this is not (yet)
         * used. In the console app (source/main.cpp), this is controlled by the
         * define WRITE_OUTPUT_TO_JSON_OBJECT. Note that not all output types
         * are implemented yet in JSONOutput, only "setup" and "eedf" sections
         * are produced at present (see source/Output.cpp).
         */
        json_type* data_out_ptr = nullptr;
        simulation.reset(new loki::Simulation(cnf, data_out_ptr));

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

