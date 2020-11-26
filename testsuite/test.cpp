#include <iostream>

#include "LoKI-B/Setup.h"
#include "LoKI-B/Simulation.h"

int main(int argc, char **argv)
{
    try
    {
        if (argc != 2)
        {
            throw std::runtime_error("Expected the input file as the single argument.");
        }
        std::string fileName(argv[1]);
        const loki::json_type cnf = loki::read_json_from_file(fileName);
        std::unique_ptr<loki::Simulation> simulation(new loki::Simulation(cnf));
//        simulation->m_obtainedResults.addListener(handleResults);
        simulation->run();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}

