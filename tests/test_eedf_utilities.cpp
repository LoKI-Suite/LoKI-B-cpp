/** \file
 *
 *  Unit tests for the code in EedfUtilities.h.
 *
 *  \author Jan van Dijk
 *  \date   19 June 2024
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/EedfUtilities.h"

#include "tests/TestUtilities.h"

using namespace loki;

void do_test(const Grid& grid, double s, double TeV, double tolerance)
{
    std::cout << "Testing s = " << s << ", TeV = " << TeV << std::endl;
    const Vector eedf_norm = makePrescribedEDF(grid,s,TeV,true);
    const Vector eedf_raw  = makePrescribedEDF(grid,s,TeV,false);
    const Vector relativeError = (eedf_norm - eedf_raw).cwiseQuotient(eedf_raw).cwiseAbs();
    test_expr( relativeError.mean() <  tolerance );
    test_expr( (getMeanEnergy(eedf_norm,grid)-1.5*TeV) < tolerance*TeV );
}

int main()
{
    using namespace loki;

    const unsigned nCells = 1000;
    const double TeV  = 1;  // eV
    const double uMax = 10; // eV

    // empirical value for the given nCells.
    const double tolerance = 1e-3;

    const Grid grid1(nCells, uMax);

    do_test(grid1,1,TeV,tolerance);
    do_test(grid1,2,TeV,tolerance);
    do_test(grid1,3,TeV,tolerance);

    // create a non-uniform grid.
    Vector faceDistribution(nCells+1);
    for (Grid::Index i=0; i!=faceDistribution.size(); ++i)
    {
        const double x = i*1.0/(faceDistribution.size()-1);
        faceDistribution[i] = (std::exp(x)-1)/(std::exp(1)-1);
    }
    const Grid grid2(faceDistribution, uMax);

    do_test(grid2,1,TeV,tolerance);
    do_test(grid2,2,TeV,tolerance);
    do_test(grid2,3,TeV,tolerance);

    test_report;
    return nerrors;
}
