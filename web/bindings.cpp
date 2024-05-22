#include <cstdint>
#include <iostream>

#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Log.h"
#include "LoKI-B/Output.h"
#include "LoKI-B/Simulation.h"
#include <chrono>
#include <exception>
#include <sstream>

#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/em_asm.h>

namespace loki
{
namespace web
{

int run(std::string file_contents, emscripten::val callback, emscripten::val output_callback)
{
    try
    {
        auto begin = std::chrono::high_resolution_clock::now();

        std::stringstream ss(file_contents);
        const json_type cnf = read_json_from_stream(ss);

        /* Create a json object that holds the output and set up a JSONOutput
         * object that will populate the JSON data object.
         */
        json_type data_out;

        std::unique_ptr<Simulation> simulation(new Simulation("", cnf));
        std::unique_ptr<Output> output(new JsonOutput(data_out, cnf, &simulation->m_workingConditions));
        /** \todo Perhaps the above should be controlled by
         * cnf.at("output").at("isOn"). \todo Now that JsonOutput works, we have
         * two ouput mechanisms in place: handleResults and handleJSONOutput. I
         * think we need only one. It may be useful to send
         * intermediate/incremental output to JS, and let that concatenate the
         * bits and pieces.
         */
        simulation->obtainedResults().addListener(
            [callback](const Grid &grid, const Vector &eedf, const WorkingConditions &wc, const Power &power,
                       const EedfCollisionDataMixture &collData, const SwarmParameters &swarmParameters,
                       const Vector *firstAnisotropy) {
                callback(wc.reducedField(),
                         reinterpret_cast<uintptr_t>(grid.getCells().data()), grid.getCells().size(),
                         reinterpret_cast<uintptr_t>(eedf.data()), eedf.size());
            });
        simulation->obtainedResults().addListener(&loki::Output::saveCycle, output.get());

        simulation->run();
        auto end = std::chrono::high_resolution_clock::now();
        Log<Message>::Notify("Simulation finished, elapsed time = ",
                             std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count(), "mus");

        // generate output
        const std::string msg = data_out.dump();
        output_callback(msg);

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}

} // namespace web
} // namespace loki

EMSCRIPTEN_BINDINGS(loki)
{
    emscripten::function("run", &loki::web::run);
}
