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
void write_grid(std::ostream& os, const loki::Vector nodeDistribution, const loki::Vector& eedf)
{
    for (loki::Grid::Index k = 0; k < nodeDistribution.size(); ++k)
    {
        os << nodeDistribution[k] << '\t' << eedf(k)  << std::endl;
    }
}

void write_eedfAnalytical(std::ostream& os, const loki::Grid& grid, const loki::Vector& eedf, double e)
{
    for (loki::Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        os << grid.getCells()(k)/e << '\t' << eedf(k)<< std::endl;
    }
}

loki::Vector analyticalSolution(loki::Grid const &grid)
{
    using namespace loki;
    double Td = 1e-21;
    double e = Constant::electronCharge;

    double uMax = 30.*e;
    int nCells = grid.nCells();
    const double eon = 10*Td;

    double u0 = 10.*e;
    double sigma0 = 1e-21;
    double Q = 1e-19;
    double mM = 2.5e-5;
    double W = e/ std::sqrt(6*mM) * eon/Q;

    Grid::Vector nodeDistribution = grid.getNodes()/grid.uMax();
    
    Grid grid1(nodeDistribution, uMax, false);
    std::ofstream ofs("grid1.dat");
    write_grid(ofs,grid1.getNodes()/e,Vector::Ones(nCells+1));
    Vector f = Vector::Zero(nCells);
    Vector fhg = Vector::Zero(nCells);
    Vector fIhg = Vector::Zero(nCells);

    Grid::Index findIndex = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), u0) - grid1.getCells().begin());
    
    for (Grid::Index k = 0; k < grid1.nCells(); ++k)
        {
            // Analytical solution
            fhg[k] = std::exp(-1*std::pow(u0/W,2)/2 * (std::pow(grid1.getCell(k)/u0,2) - std::pow(u0/u0,2)));
            
            
            if (k < findIndex)
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
    write_eedfAnalytical(ofs1,grid1,f,e);
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
    loki::Vector solution = analyticalSolution(grid);
    
    if (grid.getCells().size() != eedf.size())
    {
        throw std::runtime_error("error");
    }
    std::ofstream ofs("eedfinelastic.dat");
    write_eedf(ofs,grid,eedf);
    
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
        int n = 600;
        int n1= 200;
        
        int nCells =  n + n1;
        double U0 = 10.;
        double duPeak = 1e-1;
        double Umax = 30.;
        double deltau = 2*duPeak/n1;
        double sigma0 = 1e-21;
        std::ifstream ifs("input/JSON/NonUniform/nonuniformDelta.input.json");
        auto j = nlohmann::json::parse(ifs);
        Vector part1 = Vector::LinSpaced(n*U0/Umax, 0.0/Umax, (U0 - duPeak - 0.1*duPeak)/Umax);
        Vector part2 = Vector::LinSpaced(n1 + 1, (U0 - duPeak)/Umax, (U0 + duPeak)/Umax);
        Vector part3 = Vector::LinSpaced(n*(Umax - U0)/Umax, (U0 + duPeak + 0.1*duPeak)/Umax, Umax/Umax);

        Vector total(nCells +1);
        total << part1, part2, part3;
        std::sort(total.begin(), total.end());
        std::ofstream ofs("grid.dat");
        write_grid(ofs,total,Vector::Ones(nCells+1));
        j["electronKinetics"]["numerics"]["energyGrid"]["nonuniformGrid"]["nodeDistribution"] = total; 
        j["electronKinetics"]["numerics"]["energyGrid"]["nonuniformGrid"]["maxEnergy"] = Umax;

        j["electronKinetics"]["mixture"]["processes"][1]["threshold"] = U0 - deltau;
        j["electronKinetics"]["mixture"]["processes"][1]["data"][0][0]  = 0;
        j["electronKinetics"]["mixture"]["processes"][1]["data"][1][0]  = U0 - deltau;
        j["electronKinetics"]["mixture"]["processes"][1]["data"][2][0]  = U0;
        j["electronKinetics"]["mixture"]["processes"][1]["data"][3][0]  = U0 + deltau;
        j["electronKinetics"]["mixture"]["processes"][1]["data"][4][0]  = 1e3;
        j["electronKinetics"]["mixture"]["processes"][1]["data"][0][1]  = 0;
        j["electronKinetics"]["mixture"]["processes"][1]["data"][1][1]  = 0;
        j["electronKinetics"]["mixture"]["processes"][1]["data"][2][1]  = sigma0*U0/(deltau);
        j["electronKinetics"]["mixture"]["processes"][1]["data"][3][1]  = 0;
        j["electronKinetics"]["mixture"]["processes"][1]["data"][4][1]  = 0;

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
