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
#include "LoKI-B/Environment.h"
#include "LoKI-B/Simulation.h"
#include "tests/TestUtilities.h"
#include "source/Operators.cpp"
#include <iostream>
#include <fstream>

/// \todo can be removed once the json above is patched. See loki::legacyToJSON below.
#include "LoKI-B/LegacyToJSON.h"


int nTot = 3000;
double fraction = 0.15;

loki::Vector analyticalSolutionTwoSinglePeak(loki::Grid const &grid)
{
    using namespace loki;
    double Td = 1e-21;
    double e = Constant::electronCharge;

    double uMax = 30.*e;
    int nCells = grid.nCells();
    const double eon = 10*Td;


    double u0 = 10.*e;
    double u1 = 15.*e;
    double dupeak0 = 2*1e-1*e/(grid.nCells()*fraction);
    // double dupeak0 = uMax/nCells;
    double sigma0 = 1e-22;
    double sigma1 = 1e-19;
    double Q = 1e-19;
    double mM = 2.5e-5;
    double W = e/ std::sqrt(6*mM) * eon/Q;

    Grid::Vector nodeDistribution = grid.getNodes()/grid.uMax();
    
    Grid grid1(nodeDistribution, uMax, false);
    Vector f = Vector::Zero(nCells);
    Vector fhg = Vector::Zero(nCells);
    Vector g = Vector::Zero(nCells);

    Grid::Index findIndex1 = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), dupeak0) - grid1.getCells().begin());
    Grid::Index findIndex2 = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), u0) - grid1.getCells().begin());
    Grid::Index findIndex3 = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), u1) - grid1.getCells().begin());

    for (Grid::Index k = (grid.nCells()-1); (k >-1) ; --k)
        {
            // Analytical solution
            fhg[k] = std::exp(-1*std::pow(1/W,2)/2 * (std::pow(grid1.getCell(k),2) - std::pow(u1,2)));
            
            
            if (k < findIndex1)
            {
                g[k] = g[findIndex1];
            } else if (k < findIndex2 && k >= findIndex1)
            {
                g[k] = -1/(4*mM)*(std::pow(u0/W,2)*(sigma0/Q)*(fhg[findIndex2]*(1+g[findIndex2]))/(fhg[findIndex3]*(1 + g[findIndex3])) + std::pow(u1/W,2)*sigma1/Q)*std::exp(-std::pow(u1/W,2)/2)*(-std::expint(std::pow(u0/W,2)/2)+ std::expint(std::pow(grid1.getCell(k)/W,2)/2)) + g[findIndex2];
            } else if (k < findIndex3 && k>= findIndex2) 
            {
                g[k] = -1/(4*mM)*(sigma1/Q)*std::pow(u1/W,2)*std::exp(-1*std::pow(u1/W,2)/2)*(-std::expint(std::pow(u1/W,2)/2)+std::expint(std::pow(grid1.getCell(k)/W,2)/2));
            } else 
            {
                g[k] = 0;
            }

            f[k] = fhg[k]*(1+g[k]);
            
        }
    f /= f.dot(grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells())); 
    f *= std::pow(e,1.5);
    return f;
}

loki::Vector analyticalSolutionSinglePeak(loki::Grid const &grid)
{
    using namespace loki;
    double Td = 1e-21;
    double e = Constant::electronCharge;

    double uMax = 30.*e;
    int nCells = grid.nCells();
    const double eon = 10*Td;

    double u0 = 10.*e ;
    double sigma0 = 1e-21;
    double Q = 1e-19;
    double mM = 2.5e-5;
    double W = e/ std::sqrt(6*mM) * eon/Q;

    Grid::Vector nodeDistribution = grid.getNodes()/grid.uMax();
    
    Grid grid1(nodeDistribution, uMax, false);
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
    return f;
}

loki::Vector analyticalSolutionDoublePeak(loki::Grid const &grid)
{
    using namespace loki;
    double Td = 1e-21;
    double e = Constant::electronCharge;

    double uMax = 30.*e;
    int nCells = grid.nCells();
    const double eon = 10*Td;

    double u0 = 10.*e;
    double u1 = 15.*e;
    double dupeak0 = 2*1e-1*e/(grid.nCells()*fraction);
    double sigma0 = 1e-22;
    double sigma1 = 1e-19;
    double Q = 1e-19;
    double mM = 2.5e-5;
    double W = e/ std::sqrt(6*mM) * eon/Q;

    Grid::Vector nodeDistribution = grid.getNodes()/grid.uMax();
    
    Grid grid1(nodeDistribution, uMax, false);
    Vector f = Vector::Zero(nCells);
    Vector fhg = Vector::Zero(nCells);
    Vector g = Vector::Zero(nCells);

    Grid::Index findIndex1 = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), dupeak0) - grid1.getCells().begin());
    Grid::Index findIndex2 = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), (u1-u0)) - grid1.getCells().begin());
    Grid::Index findIndex3 = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), u0) - grid1.getCells().begin());
    Grid::Index findIndex4 = (std::upper_bound(grid1.getCells().begin(),grid1.getCells().end(), u1) - grid1.getCells().begin());

    for (Grid::Index k = (grid.nCells()); (k >-1) ; --k)
        {
            // Analytical solution
            fhg[k] = std::exp(-1.0*std::pow(1.0/W,2)/2.0 * (std::pow(grid1.getCell(k),2) - std::pow(u1,2)));
            
            
            if (k < findIndex1)
            {
                g[k] = g[findIndex1];
            } else if (k < findIndex2 && k >= findIndex1)
            {
                g[k] = 1./(4.*mM)*(sigma0/Q)*std::pow(u0/W,2)*(fhg[findIndex3]*(1+g[findIndex3]))/(fhg[findIndex4]*(1+g[findIndex4]))* std::exp(-1*std::pow(u1/W,2)/2) * (std::expint(std::pow((u1-u0)/W,2)/2.) - std::expint(std::pow(grid1.getCell(k)/W,2)/2.)) + g[findIndex2];
            } else if (k < findIndex3 && k>= findIndex2) 
            {
                g[k] = 1./(4.*mM)*(std::pow(u0/W,2)*(sigma0/Q)*(fhg[findIndex3]*(1+g[findIndex3]))/(fhg[findIndex4]*(1+g[findIndex4])) + std::pow(u1/W,2)*(sigma1/Q))*std::exp(-1*std::pow(u1/W,2)/2.)*(std::expint(std::pow(u0/W,2)/2.) - std::expint(std::pow(grid1.getCell(k)/W,2)/2.)) + g[findIndex3];
            } else if (k < findIndex4 && k >= findIndex3)
            {
                g[k] = 1./(4.*mM)*(sigma1/Q)*std::pow(u1/W,2)* std::exp(-1*std::pow(u1/W,2)/2) *(std::expint(std::pow(u1/W,2)/2) - std::expint(std::pow(grid1.getCell(k)/W,2)/2));
            } else 
            {
                g[k] = 0;
            }

            f[k] = fhg[k]*(1+g[k]);
            
        }
    f /= f.dot(grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells())); 
    f *= std::pow(e,1.5);
    return f;
}

void checkSinglePeak(const loki::Grid &grid,
    const loki::Vector &eedf,
    const loki::WorkingConditions &wc,
    const loki::Power &power,
    const loki::EedfCollisionDataMixture& gases,
    const loki::SwarmParameters &swarmParameters,
    const loki::Vector *firstAnisotropy
)
{
    using namespace loki;
    loki::Vector solution = analyticalSolutionSinglePeak(grid);
    
    if (grid.getCells().size() != eedf.size())
    {
        throw std::runtime_error("error");
    }
    
    double meanRelativeError = ((solution - eedf).cwiseQuotient(solution)).cwiseAbs().mean();
    test_expr(meanRelativeError < 0.50);
}

void checkTwoSinglePeak(const loki::Grid &grid,
    const loki::Vector &eedf,
    const loki::WorkingConditions &wc,
    const loki::Power &power,
    const loki::EedfCollisionDataMixture& gases,
    const loki::SwarmParameters &swarmParameters,
    const loki::Vector *firstAnisotropy
)
{
    using namespace loki;
    loki::Vector solution = analyticalSolutionTwoSinglePeak(grid);
    if (grid.getCells().size() != eedf.size())
    {
        throw std::runtime_error("error");
    }
    
    double meanRelativeError = ((solution - eedf).cwiseQuotient(solution)).cwiseAbs().mean();
    test_expr(meanRelativeError < 0.50);
}

void checkDoublePeak(const loki::Grid &grid,
    const loki::Vector &eedf,
    const loki::WorkingConditions &wc,
    const loki::Power &power,
    const loki::EedfCollisionDataMixture& gases,
    const loki::SwarmParameters &swarmParameters,
    const loki::Vector *firstAnisotropy
)
{
    using namespace loki;
    loki::Vector solution = analyticalSolutionDoublePeak(grid);
    
    if (grid.getCells().size() != eedf.size())
    {
        throw std::runtime_error("error");
    }
    
    double meanRelativeError = ((solution - eedf).cwiseQuotient(solution)).cwiseAbs().mean();
    test_expr(meanRelativeError < 0.50);
}

nlohmann::json twoSingleDeltaPeaks(int nCells, double frac)
{
    using namespace loki;
    
    int n1 = nCells * frac;
    int n =  nCells - 2*n1;
    double U0 = 10.;
    double U1 = 15.;
    double duPeak = 1e-1;
    double Umax = 30.;
    double deltau = 2*duPeak/n1;

    // If grid is uniform : 
    // double deltau = Umax/nCells;

    double sigma0 = 1e-22;
    double sigma1 = 1e-19;
    const std::string fname = getEnvironmentVariable("LOKI_TEST_INPUT_DIR") + "/../input/JSON/Delta/nonuniform-two-single.json";
    std::ifstream ifs(fname);
    if (!ifs)
    {
        throw std::runtime_error("Could not open file '" + fname + "' for reading.");
    }
    auto j = nlohmann::json::parse(ifs);
    /// \todo the json literal above is still 'old JSON'. patch it for the time being.
    j = loki::legacyToJSON(j);
    double zero = 0.03 * deltau;
    Vector part1 = Vector::LinSpaced(n*U0/Umax, 0.0/Umax, (U0 - duPeak - 0.1*duPeak)/Umax);
    Vector part2 = Vector::LinSpaced(n1, (U0 - duPeak)/Umax, (U0 + duPeak)/Umax);
    Vector part3 = Vector::LinSpaced(n*(U1 - U0)/Umax, (U0 + duPeak + 0.1*duPeak)/Umax, (U1 - duPeak - 0.1*duPeak)/Umax);
    Vector part4 = Vector::LinSpaced(n1, (U1 - duPeak)/Umax, (U1 + duPeak)/Umax);
    Vector part5 = Vector::LinSpaced(n*(Umax - U1)/Umax, (U1 + duPeak + 0.1*duPeak)/Umax, Umax/Umax);

    Vector nodeDistribution(nCells +1);
    nodeDistribution << zero, part1, part2, part3, part4, part5;
    std::sort(nodeDistribution.begin(), nodeDistribution.end());
    j["electronKinetics"]["numerics"]["energyGrid"]["nonuniformGrid"]["nodeDistribution"] = nodeDistribution; 
    j["electronKinetics"]["numerics"]["energyGrid"]["nonuniformGrid"]["maxEnergy"] = Umax;

    // Uniform grid : 
    // j["electronKinetics"]["numerics"]["energyGrid"]["maxEnergy"] = Umax;
    // j["electronKinetics"]["numerics"]["energyGrid"]["cellNumber"] = nCells;

    auto &firstProcessInfo = j["electronKinetics"]["mixture"]["processes"][1]["info"][0];

    firstProcessInfo["threshold"] = U0 - deltau;
    firstProcessInfo["data"]["values"][0][0]  = 0;
    firstProcessInfo["data"]["values"][1][0]  = U0 - deltau;
    firstProcessInfo["data"]["values"][2][0]  = U0;
    firstProcessInfo["data"]["values"][3][0]  = U0 + deltau;
    firstProcessInfo["data"]["values"][4][0]  = 1e3;
    firstProcessInfo["data"]["values"][0][1]  = 0;
    firstProcessInfo["data"]["values"][1][1]  = 0;
    firstProcessInfo["data"]["values"][2][1]  = sigma0*U0/(deltau);
    firstProcessInfo["data"]["values"][3][1]  = 0;
    firstProcessInfo["data"]["values"][4][1]  = 0;

    auto &secondProcessInfo = j["electronKinetics"]["mixture"]["processes"][2]["info"][0];

    secondProcessInfo["threshold"] = U1 - deltau;
    secondProcessInfo["data"]["values"][0][0]  = 0;
    secondProcessInfo["data"]["values"][1][0]  = U1 - deltau;
    secondProcessInfo["data"]["values"][2][0]  = U1;
    secondProcessInfo["data"]["values"][3][0]  = U1 + deltau;
    secondProcessInfo["data"]["values"][4][0]  = 1e3;
    secondProcessInfo["data"]["values"][0][1]  = 0;
    secondProcessInfo["data"]["values"][1][1]  = 0;
    secondProcessInfo["data"]["values"][2][1]  = sigma1*U1/(deltau);
    secondProcessInfo["data"]["values"][3][1]  = 0;
    secondProcessInfo["data"]["values"][4][1]  = 0;

    return j;
}

nlohmann::json doubleDeltaPeaks(int nCells, double frac)
{
    using namespace loki;

    int n1 = nCells * frac;
    int n =  nCells - 2*n1;
    double U0 = 10.;
    double alpha = 1.5;
    double U1 = U0*alpha;
    double duPeak1 = 1e-1;
    double duPeak2 = 1e-1;
    double Umax = 30.;

    double deltau1 = 2*duPeak1/n1;
    double deltau2 = 2*duPeak2/n1;

    // If grid is uniform : 
    // double deltau1 = Umax/nCells;
    // double deltau2 = Umax/nCells;
    
    double sigma0 = 1e-22;
    double sigma1 = 1e-19;

    const std::string fname = getEnvironmentVariable("LOKI_TEST_INPUT_DIR") + "/../input/JSON/Delta/nonuniform-double.json";
    std::ifstream ifs(fname);
    if (!ifs)
    {
        throw std::runtime_error("Could not open file '" + fname + "' for reading.");
    }
    auto j = nlohmann::json::parse(ifs);
    /// \todo the json literal above is still 'old JSON'. patch it for the time being.
    j = loki::legacyToJSON(j);
    Vector part1 = Vector::LinSpaced(n*U0/Umax, 0.0/Umax, (U0 - duPeak1 - 0.1*duPeak1)/Umax);
    Vector part2 = Vector::LinSpaced(n1, (U0 - duPeak1)/Umax, (U0 + duPeak1)/Umax);
    Vector part3 = Vector::LinSpaced(n*(U1 - U0)/Umax + 1, (U0 + duPeak1 + 0.1*duPeak1)/Umax, (U1 - duPeak2 - 0.1*duPeak2)/Umax);
    Vector part4 = Vector::LinSpaced(n1, (U1 - duPeak2)/Umax, (U1 + duPeak2)/Umax);
    Vector part5 = Vector::LinSpaced(n*(Umax - U1)/Umax, (U1 + duPeak2 + 0.1*duPeak2)/Umax, Umax/Umax);


    Vector nodeDistribution(nCells +1);
    nodeDistribution << part1, part2, part3, part4, part5;
    std::sort(nodeDistribution.begin(), nodeDistribution.end());
    j["electronKinetics"]["numerics"]["energyGrid"]["nonuniformGrid"]["nodeDistribution"] = nodeDistribution; 
    j["electronKinetics"]["numerics"]["energyGrid"]["nonuniformGrid"]["maxEnergy"] = Umax;

    // Uniform grid : 
    // j["electronKinetics"]["numerics"]["energyGrid"]["maxEnergy"] = Umax;
    // j["electronKinetics"]["numerics"]["energyGrid"]["cellNumber"] = nCells;

    auto &processInfo = j["electronKinetics"]["mixture"]["processes"][1]["info"][0];

    processInfo["threshold"] = U0 - deltau1;
    processInfo["data"]["values"][0][0]  = 0;
    processInfo["data"]["values"][1][0]  = U0 - deltau1;
    processInfo["data"]["values"][2][0]  = U0;
    processInfo["data"]["values"][3][0]  = U0 + deltau1;
    processInfo["data"]["values"][4][0]  = U1 - deltau2;
    processInfo["data"]["values"][5][0]  = U1;
    processInfo["data"]["values"][6][0]  = U1 + deltau2;
    processInfo["data"]["values"][7][0]  = 1e3;

    processInfo["data"]["values"][0][1]  = 0;
    processInfo["data"]["values"][1][1]  = 0;
    processInfo["data"]["values"][2][1]  = sigma0*U0/(deltau1);
    processInfo["data"]["values"][3][1]  = 0;
    processInfo["data"]["values"][4][1]  = 0;
    processInfo["data"]["values"][5][1]  = sigma1*U1/(deltau2);
    processInfo["data"]["values"][6][1]  = 0;
    processInfo["data"]["values"][7][1]  = 0;

    return j;
}

nlohmann::json singleDeltaPeak(int nCells, double frac)
{
    using namespace loki;
    double duPeak = 1e-1;
    int n1 = nCells * frac;
    int n =  nCells - n1;
    
    double U0 = 10.;
    double Umax = 30.;
    double deltau = 2*duPeak/(n1);

    // If grid is uniform
    // double deltau = Umax/nCells;

    double sigma0 = 1e-21;
    const std::string fname = getEnvironmentVariable("LOKI_TEST_INPUT_DIR") + "/../input/JSON/Delta/nonuniform.json";
    std::ifstream ifs(fname);
    if (!ifs)
    {
        throw std::runtime_error("Could not open file '" + fname + "' for reading.");
    }
    auto j = nlohmann::json::parse(ifs);
    /// \todo the json literal above is still 'old JSON'. patch it for the time being.
    j = loki::legacyToJSON(j);
    double zero = deltau * 0.03;
    Vector part1 = Vector::LinSpaced(n*U0/Umax, 0.0/Umax, (U0 - duPeak - 0.1*duPeak)/Umax);
    Vector part2 = Vector::LinSpaced(n1, (U0 - duPeak)/Umax, (U0 + duPeak)/Umax);
    Vector part3 = Vector::LinSpaced(n*(Umax - U0)/Umax, (U0 + duPeak + 0.1*duPeak)/Umax, Umax/Umax);

    Vector nodeDistribution(nCells +1);
    nodeDistribution << zero, part1, part2, part3;
    std::sort(nodeDistribution.begin(), nodeDistribution.end());
    j["electronKinetics"]["numerics"]["energyGrid"]["nonuniformGrid"]["nodeDistribution"] = nodeDistribution; 
    j["electronKinetics"]["numerics"]["energyGrid"]["nonuniformGrid"]["maxEnergy"] = Umax;
    
    // Uniform grid : 
    // j["electronKinetics"]["numerics"]["energyGrid"]["maxEnergy"] = Umax;
    // j["electronKinetics"]["numerics"]["energyGrid"]["cellNumber"] = nCells;

    auto &processInfo = j["electronKinetics"]["mixture"]["processes"][1]["info"][0];

    processInfo["threshold"] = U0- deltau;
    processInfo["data"]["values"][0][0]  = 0;
    processInfo["data"]["values"][1][0]  = U0 - deltau;
    processInfo["data"]["values"][2][0]  = U0;
    processInfo["data"]["values"][3][0]  = U0 + deltau;
    processInfo["data"]["values"][4][0]  = 1e3;
    processInfo["data"]["values"][0][1]  = 0;
    processInfo["data"]["values"][1][1]  = 0;
    processInfo["data"]["values"][2][1]  = sigma0*U0/(deltau);
    processInfo["data"]["values"][3][1]  = 0;
    processInfo["data"]["values"][4][1]  = 0;

    return j;
}

int main()
{
    using namespace loki;
    try
    {
        auto i = singleDeltaPeak(nTot, fraction);
        std::unique_ptr<loki::Simulation> simulationSingle(new loki::Simulation("", i));
        simulationSingle->obtainedResults().addListener(checkSinglePeak);
        simulationSingle->run();

        auto j = twoSingleDeltaPeaks(nTot, fraction);
        std::unique_ptr<loki::Simulation> simulationTwoSingle(new loki::Simulation("", j));
        simulationTwoSingle->obtainedResults().addListener(checkTwoSinglePeak);
        simulationTwoSingle->run();

        auto m = doubleDeltaPeaks(nTot, fraction);
        // std::unique_ptr<loki::Simulation> simulationDouble(new loki::Simulation(m));
        // simulationDouble->obtainedResults().addListener(checkDoublePeak);
        // simulationDouble->run();
    }
    catch (const std::exception &exc)
    {
        std::cerr << exc.what() << std::endl;
        return 1;
    }

    test_report;
    return nerrors;
}
