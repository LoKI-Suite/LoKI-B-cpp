/** \file
 *
 *  Unit tests for Grid.h
 *
 *  \author Jan van Dijk
 *  \date   September 2023
 */

#include "LoKI-B/Grid.h"

#include "tests/TestUtilities.h"

int main()
{
    using namespace loki;
    const double eps = std::numeric_limits<double>::epsilon();

    Grid grid1(5,5);
    test_expr( std::abs(grid1.duCells().sum()-grid1.uMax()) < eps*grid1.uMax() );
    test_expr( std::abs(grid1.duNodes().sum()-grid1.uMax()) < eps*grid1.uMax() );
    grid1.updateMaxEnergy(10);
    test_expr( std::abs(grid1.duCells().sum()-grid1.uMax()) < eps*grid1.uMax() );
    test_expr( std::abs(grid1.duNodes().sum()-grid1.uMax()) < eps*grid1.uMax() );

    Grid::Vector ls(6); ls << 0.0,0.2,0.4,0.6,0.8,1.0;
    Grid grid2(ls,10,true); // true: we promise that ls is indeed equidistant
    test_expr( std::abs(grid2.duCells().sum()-grid2.uMax()) < eps*grid2.uMax() );
    test_expr( std::abs(grid2.duNodes().sum()-grid2.uMax()) < eps*grid2.uMax() );

    test_expr( (grid1.getNodes()-grid2.getNodes()).squaredNorm()<eps*grid1.uMax()/grid1.getNodes().size() );

    test_report;
    return nerrors;
}
