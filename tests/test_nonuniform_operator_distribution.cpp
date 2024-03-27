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

void test_distribution(loki::Vector nodes, unsigned int sourceIdx, bool reverse = false)
{
    int rev = reverse ? 1 : -1;
    using namespace loki;
    int uMax = 100;
    double threshold = 50.;
    const double eps = std::numeric_limits<double>::epsilon();

    Grid grid(nodes, uMax, false);
    std::vector<std::tuple<int, double>> alpha;
    alpha = getOperatorDistribution(grid, threshold, grid.getCell(sourceIdx), sourceIdx, reverse);

    double alphaTot = 0;
    double energyTarget = 0;
    for (int i = 0; i < int(alpha.size()); i++)
    {
        alphaTot += std::get<1>(alpha[i]);
        energyTarget += grid.getCell(std::get<0>(alpha[i])) * std::get<1>(alpha[i]);
    }
    test_expr(std::abs(alphaTot - 1.) <= eps);
    test_expr(std::abs(energyTarget - (grid.getCell(sourceIdx) + rev*threshold)) <= eps);
}

int main()
{
    using namespace loki;

    Vector nodeDistribution1(8);
    nodeDistribution1 << 0.0, 0.1, 0.2, 0.3, 0.5, 0.65, 0.75, 1.0;
    test_distribution(nodeDistribution1, 5);

    Vector nodeDistribution2(7);
    nodeDistribution2 << 0.0, 0.1, 0.3, 0.5, 0.65, 0.78, 1.0;
    test_distribution(nodeDistribution2, 4);

    Vector nodeDistribution3(12);
    nodeDistribution3 << 0.0, 0.1, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40, 0.5, 0.63, 0.87, 1.0;
    test_distribution(nodeDistribution3, 9);

    Vector nodeDistribution4(7);
    nodeDistribution4 << 0.0, 0.12, 0.25, 0.5, 0.6, 0.8, 1.0;
    test_distribution(nodeDistribution4, 1, true);

    Vector nodeDistribution5(8);
    nodeDistribution5 << 0.0, 0.15, 0.28, 0.5, 0.6, 0.7, 0.8, 1.0;
    test_distribution(nodeDistribution5, 1, true);

    Vector nodeDistribution6(12);
    nodeDistribution6 << 0.0, 0.13, 0.38, 0.5, 0.6, 0.65, 0.70, 0.75, 0.8, 0.85, 0.90, 1.0;
    test_distribution(nodeDistribution6, 1, true);

    test_report;
    return nerrors;
}