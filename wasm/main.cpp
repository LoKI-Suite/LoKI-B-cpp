#include <iostream>

#include <Setup.h>
#include <Simulation.h>
#include <LinearAlgebra.h>
#include <chrono>
#include <Log.h>

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::WorkingConditions &wc,
                   const loki::Power &power, const std::vector<loki::EedfGas *> &gasses,
                   const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients,
                   const loki::Vector &firstAnisotropy);

void handleExistingOutputPath(std::string &folder);

void plot(const std::string &title, const std::string &xlabel, const std::string &ylabel,
          const loki::Vector &x, const loki::Vector &y);

int main()
{
    std::cout << "Loaded loki-b module." << std::endl;
    return 0;
}
extern "C"
{
    int run(char *a)
    {
        std::string fileName = a;

        try
        {
            loki::Setup setup;

            if (!setup.parseFile(fileName))
            {
                return 1;
            }

            loki::Simulation simulation(setup);

            simulation.obtainedResults.addListener(handleResults);
            simulation.outputPathExists.addListener(handleExistingOutputPath);

            simulation.run();
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }

        return 0;
    }
}

void handleResults(const loki::Grid &grid, const loki::Vector &eedf, const loki::WorkingConditions &wc,
                   const loki::Power &power, const std::vector<loki::EedfGas *> &gasses,
                   const loki::SwarmParameters &swarmParameters,
                   const std::vector<loki::RateCoefficient> &rateCoefficients,
                   const std::vector<loki::RateCoefficient> &extraRateCoefficients,
                   const loki::Vector &firstAnisotropy)
{
    plot("Eedf", "Energy (eV)", "Eedf (Au)", grid.getCells(), eedf);
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
