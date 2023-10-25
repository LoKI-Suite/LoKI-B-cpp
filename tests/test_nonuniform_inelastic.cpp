/** \file
 *
 *  Unit tests for nonuniform inelastic operator
 *
 *  \author Jop Hendrikx
 *  \date   October 2023
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/EedfUtilities.h"
#include "LoKI-B/EedfCollisions.h"
#include "LoKI-B/Simulation.h"
#include "tests/TestUtilities.h"
#include "source/Operators.cpp"
#include <iostream>
#include <fstream>

void write_eedf(std::ostream& os, const loki::Grid& grid, const loki::Vector& eedf)
{
    for (loki::Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        os << grid.getCells()(k) << '\t' << eedf(k)  << std::endl;
    }
}

void write_eedfAnalytical(std::ostream& os, const loki::Grid& grid, const loki::Vector& eedf, double e)
{
    for (loki::Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        os << grid.getCells()(k)/e << '\t' << eedf(k)<< std::endl;
    }
}

loki::Vector analyticalSolution()
{
    using namespace loki;
    double Td = 1e-21;
    double e = Constant::electronCharge;

    double uMax = 30.*e;
    int nCells = 6000;
    const double eon = 10*Td;

    double u0 = 10.*e;
    double sigma0 = 1e-21;
    double Q = 1e-19;
    double mM = 2.5e-5;
    double W = e/ std::sqrt(6*mM) * eon/Q;

    Grid::Vector nodeDistribution = Vector::LinSpaced(nCells + 1, 0.0, 1.0);
    Grid grid1(nodeDistribution, uMax, false);
    Grid grid2(nCells, uMax);
    Vector f = Vector::Zero(nCells);
    Vector fhg = Vector::Zero(nCells);
    Vector fIhg = Vector::Zero(nCells);

    Grid::Index findIndex = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), u0) - grid1.getCells().begin());
    
    for (Grid::Index k = 0; k < grid1.nCells(); ++k)
        {
            // Analytical solution
            fhg[k] = std::exp(-1*std::pow(u0/W,2)/2 * (std::pow(grid1.getCell(k)/u0,2) - std::pow(u0/u0,2)));
            
            
            if (k <= findIndex)
            {
                fIhg[k] = 1/(mM * 4)*sigma0/Q * std::pow(u0/W, 2)*std::exp(-std::pow(u0/W,2)/2)*(-std::expint(std::pow(grid1.getCell(k)/W,2) / 2) + std::expint(std::pow(u0/W,2) / 2));
                f[k] = fhg[k]*(1+fIhg[k]);

            } else
            {
                f[k] = fhg[k];
            }
            
        }
    f /= f.dot(grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells())); 
    f *= std::pow(e,1.5); 
    std::ofstream ofs1("eedfinelasticAnalytical.dat");
    write_eedfAnalytical(ofs1,grid2,f,e);
    return f;
}

void checkRMSE(const loki::Grid &grid,
    const loki::Vector &eedf,
    const loki::WorkingConditions &wc,
    const loki::Power &power,
    const loki::EedfCollisionDataMixture& gases,
    const loki::SwarmParameters &swarmParameters,
    const loki::Vector *firstAnisotropy
)
{
    using namespace loki;
    if (grid.getCells().size() != eedf.size())
    {
        throw std::runtime_error("error");
    }
    std::ofstream ofs("eedfinelastic.dat");
    write_eedf(ofs,grid,eedf);

    Vector solution = analyticalSolution();
    Vector relativeError = ((solution - eedf).cwiseQuotient(solution).cwiseAbs());
    std::cerr << relativeError.mean() << std::endl;
    std::ofstream ofs2("eedferror.dat");
    write_eedf(ofs2,grid,relativeError);
}


int main()
{
    using namespace loki;

    try
    {
        std::ifstream ifs("input/JSON/NonUniform/delta.input.json");
        auto j = nlohmann::json::parse(ifs);

        std::unique_ptr<loki::Simulation> simulation(new loki::Simulation(j));
        simulation->m_obtainedResults.addListener(checkRMSE);
        simulation->run();
    }
    catch (const std::exception &exc)
    {
        std::cerr << exc.what() << std::endl;
        return 1;
    }

    test_report;
    return nerrors;
}
