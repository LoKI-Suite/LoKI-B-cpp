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

int run(std::string file_contents)
{
    try
    {
        auto begin = std::chrono::high_resolution_clock::now();

        std::stringstream ss(file_contents);
        const loki::json_type cnf = loki::read_json_from_stream(ss);

        std::unique_ptr<loki::Simulation> simulation(new loki::Simulation(cnf));
        /* Create a json object that holds the output and set up a JSONOutput
         * object that will populate the JSON data object.
         */
        json_type data_out;
        simulation->configureOutput(new loki::JsonOutput(data_out, cnf,
                            &simulation->m_workingConditions, &simulation->m_jobManager));
        /** \todo Perhaps the above should be controlled by cnf.at("output").at("isOn").
         *  \todo Now that JsonOutput works, we have two ouput mechanisms in place: handleResults
         *        and handleJSONOutput. I think we need only one. It may be useful to
         *        send intermediate/incremental output to JS, and let that concatenate
         *        the bits and pieces.
         */
        simulation->m_obtainedResults.addListener(handleResults);

        simulation->run();
        auto end = std::chrono::high_resolution_clock::now();
        Log<Message>::Notify("Simulation finished, elapsed time = ",
                             std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count(), "mus");

        // generate output
        const std::string msg=data_out.dump();
        EM_ASM({ handleJsonOutput($0, $1); }, msg.data(), msg.size());

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

