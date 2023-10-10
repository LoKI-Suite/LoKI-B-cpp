/** \file
 *
 *  Unit tests for Maxwell distribution
 *
 *  \author Jop Hendrikx
 *  \date   September 2023
 */

#include "LoKI-B/Grid.h"
#include "source/Operators.cpp"
#include "source/LinearAlgebra.cpp"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Gnuplot.h"
#include "LoKI-B/EedfUtilities.h"
#include <cmath>

#include "tests/TestUtilities.h"

int main()
{
    using namespace loki;

    const unsigned nCells = 1000;
    const double T = 300; 
    const double uMax = 2; // eV
    double kT = Constant::kBeV*T;

    Vector elasticCrossSection = Vector::Ones(nCells+1);

    auto faceDistribution = Grid::Vector::LinSpaced(nCells + 1, 0.0, 1.0); 
    Grid grid(faceDistribution, uMax, false); 

    ElasticOperator elasticOperator;
    SparseMatrix M(nCells, nCells);
    elasticOperator.evaluate(grid, elasticCrossSection, T, M);

    Vector eedf = Vector::Zero(nCells);
    eedf[0] = 1.;

    Matrix boltzmann = M.toDense();
    boltzmann.row(0) = grid.getCells().cwiseSqrt().cwiseProduct(grid.duCells());

    LinAlg::hessenberg(boltzmann.data(), eedf.data(), grid.nCells());
    
    // Analytical solution
    Vector eedfMaxwell = Vector::Zero(nCells);
    eedfMaxwell << makePrescribedEDF(grid,1,kT);


    // Calculate relative error
    Vector relativeError = (eedf - eedfMaxwell).cwiseQuotient(eedfMaxwell).cwiseAbs();
    test_expr((relativeError.mean() <  0.02));

    test_report;
    return nerrors;
}