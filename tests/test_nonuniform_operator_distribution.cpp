/** \file
 *
 *  Unit tests for nonuniform inelastic operator
 *
 *  \author Jop Hendrikx
 *  \date   October 2023
 */

#include "LoKI-B/Grid.h"
#include "tests/TestUtilities.h"
#include "source/Operators.cpp"

int main()
{
    using namespace loki;
    
    int uMax = 100;
    double threshold = 50.;
    const double eps = std::numeric_limits<double>::epsilon();

    int nCells1 = 7;
    int sourceIdx1 = 5;
    Vector nodeDistribution1(nCells1+1);
    nodeDistribution1 << 0.0, 0.1, 0.2, 0.3, 0.5, 0.65, 0.75, 1.0;
    Grid grid1(nodeDistribution1, uMax, false);
    std::vector<std::tuple<int, double>> alpha1;
    alpha1 = getOperatorDistribution(grid1, threshold, grid1.getCell(sourceIdx1), sourceIdx1);

    double alphaTot1 = 0;
    double energyTarget1 = 0;
    for (int i = 0; i < int(alpha1.size()); i++)
    {
        alphaTot1 += std::get<1>(alpha1[i]);
        energyTarget1 += grid1.getCell(std::get<0>(alpha1[i])) * std::get<1>(alpha1[i]);
    }
    test_expr(std::abs(alphaTot1 - 1.) <= eps);
    test_expr(std::abs(energyTarget1 - (grid1.getCell(sourceIdx1) - threshold)) <= eps);


    int nCells2 = 6;
    int sourceIdx2 = 4;
    Vector nodeDistribution2(nCells2+1);
    nodeDistribution2 << 0.0, 0.1, 0.3, 0.5, 0.65, 0.78, 1.0;
    Grid grid2(nodeDistribution2, uMax, false);
    std::vector<std::tuple<int, double>> alpha2;
    alpha2 = getOperatorDistribution(grid2, threshold, grid2.getCell(sourceIdx2), sourceIdx2);

    double alphaTot2 = 0;
    double energyTarget2 = 0;
    for (int i = 0; i < int(alpha2.size()); i++)
    {
        alphaTot2 += std::get<1>(alpha2[i]);
        energyTarget2 += grid2.getCell(std::get<0>(alpha2[i])) * std::get<1>(alpha2[i]);
    }
    test_expr(std::abs(alphaTot2 - 1.) <= eps);
    test_expr(std::abs(energyTarget2 - (grid2.getCell(sourceIdx2) - threshold)) <= eps);


    int nCells3 = 11;
    int sourceIdx3 = 9;
    Vector nodeDistribution3(nCells3+1);
    nodeDistribution3 << 0.0, 0.1, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40, 0.5, 0.63, 0.87, 1.0;
    Grid grid3(nodeDistribution3, uMax, false);
    std::vector<std::tuple<int, double>> alpha3;
    alpha3 = getOperatorDistribution(grid3, threshold, grid3.getCell(sourceIdx3), sourceIdx3);
    double alphaTot3 = 0;
    double energyTarget3 = 0;
    for (int i = 0; i < int(alpha3.size()); i++)
    {
        alphaTot3 += std::get<1>(alpha3[i]);
        energyTarget3 += grid3.getCell(std::get<0>(alpha3[i])) * std::get<1>(alpha3[i]);
    }
    test_expr(std::abs(alphaTot3 - 1.) <= eps);
    test_expr(std::abs(energyTarget3 - (grid3.getCell(sourceIdx3) - threshold)) <= eps);

    int nCells4 = 6;
    int sourceIdx4 = 1;
    Vector nodeDistribution4(nCells4+1);
    nodeDistribution4 << 0.0, 0.12, 0.25, 0.5, 0.6, 0.8, 1.0;
    Grid grid4(nodeDistribution4, uMax, false);
    std::vector<std::tuple<int, double>> alpha4;
    alpha4 = getOperatorDistribution(grid4, threshold, grid4.getCell(sourceIdx4), sourceIdx4, true);
    double alphaTot4 = 0;
    double energyTarget4 = 0;
    for (int i = 0; i < int(alpha4.size()); i++)
    {
        alphaTot4 += std::get<1>(alpha4[i]);
        energyTarget4 += grid4.getCell(std::get<0>(alpha4[i])) * std::get<1>(alpha4[i]);
    }
    test_expr(std::abs(alphaTot4 - 1.) <= eps);
    test_expr(std::abs(energyTarget4 - (grid4.getCell(sourceIdx4) + threshold)) <= eps);

    int nCells5 = 7;
    int sourceIdx5 = 1;
    Vector nodeDistribution5(nCells5+1);
    nodeDistribution5 << 0.0, 0.15, 0.28, 0.5, 0.6, 0.7, 0.8, 1.0;
    Grid grid5(nodeDistribution5, uMax, false);
    std::vector<std::tuple<int, double>> alpha5;
    alpha5 = getOperatorDistribution(grid5, threshold, grid5.getCell(sourceIdx5), sourceIdx5, true);
    double alphaTot5 = 0;
    double energyTarget5 = 0;
    for (int i = 0; i < int(alpha5.size()); i++)
    {
        alphaTot5 += std::get<1>(alpha5[i]);
        energyTarget5+= grid5.getCell(std::get<0>(alpha5[i])) * std::get<1>(alpha5[i]);
    }
    test_expr(std::abs(alphaTot5 - 1.) <= eps);
    test_expr(std::abs(energyTarget5 - (grid5.getCell(sourceIdx5) + threshold)) <= eps);

    int nCells6 = 11;
    int sourceIdx6 = 1;
    Vector nodeDistribution6(nCells6+1);
    nodeDistribution6 << 0.0, 0.13, 0.38, 0.5, 0.6, 0.65, 0.70, 0.75, 0.8, 0.85, 0.90, 1.0;
    Grid grid6(nodeDistribution6, uMax, false);
    std::vector<std::tuple<int, double>> alpha6;
    alpha6 = getOperatorDistribution(grid6, threshold, grid6.getCell(sourceIdx6), sourceIdx6, true);
    double alphaTot6 = 0;
    double energyTarget6 = 0;
    for (int i = 0; i < int(alpha6.size()); i++)
    {
        alphaTot6 += std::get<1>(alpha6[i]);
        energyTarget6+= grid6.getCell(std::get<0>(alpha6[i])) * std::get<1>(alpha6[i]);
    }
    test_expr(std::abs(alphaTot6 - 1.) <= eps);
    test_expr(std::abs(energyTarget6 - (grid6.getCell(sourceIdx6) + threshold)) <= eps);

    test_report;
    return nerrors;
}