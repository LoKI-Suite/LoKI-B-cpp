#include <iostream>

#include "LoKI-B/Setup.h"
#include "LoKI-B/Simulation.h"

void printResults(const loki::Grid &grid,
    const loki::Vector &eedf,
    const loki::WorkingConditions &wc,
    const loki::Power &power,
    const std::vector<loki::EedfGas *> &gases,
    const loki::SwarmParameters &swarmParameters,
    const std::vector<loki::RateCoefficient> &rateCoefficients,
    const std::vector<loki::RateCoefficient> &extraRateCoefficients,
    const loki::Vector &firstAnisotropy)
{
    using namespace loki;
    if (grid.getCells().size() != eedf.size())
    {
        throw std::runtime_error("error");
    }
    for (unsigned i = 0; i < grid.getCells().size(); i++)
    {
        std::cout << grid.getCells()[i] << "\t" << eedf[i] << std::endl;
    }
}

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
        simulation->m_obtainedResults.addListener(printResults);
        simulation->run();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}

