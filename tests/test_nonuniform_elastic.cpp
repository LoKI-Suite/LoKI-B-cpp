/** \file
 *
 *  Unit tests for nonuniform elastic operator
 *
 *  \author Jop Hendrikx
 *  \date   September 2023
 */

#include "LoKI-B/Grid.h"
#include "source/Operators.cpp"
#include "LoKI-B/Constant.h"

#include "tests/TestUtilities.h"

int main()
{
    using namespace loki;

    const unsigned nCells = 5;
    const double uMax = 10; // eV
    const double T = 300;
    Vector elasticCrossSection = Vector::Ones(nCells+1);

    auto faceDistribution = Grid::Vector::LinSpaced(nCells + 1, 0.0, 1.0); 
    Grid grid1(faceDistribution, uMax, false); 
    Grid grid2(nCells, uMax);

    ElasticOperator elasticOperator;
    SparseMatrix M1(nCells, nCells);
    SparseMatrix M2(nCells, nCells);
    
    elasticOperator.evaluate(grid1, elasticCrossSection, T, M1);
    elasticOperator.evaluate(grid2, elasticCrossSection, T, M2);
    test_expr(M1.isApprox(M2));

    Vector eedf1 = Vector::Zero(nCells);
    eedf1[0] = 1.;

    Vector eedf2 = Vector::Zero(nCells);
    eedf2[0] = 1.;

    Matrix boltzmann1 = M1.toDense();
    boltzmann1.row(0) = grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells());

    Matrix boltzmann2 = M2.toDense();
    boltzmann2.row(0) = grid2.getCells().cwiseSqrt()*grid2.du();

    LinAlg::hessenberg(boltzmann1.data(), eedf1.data(), grid1.nCells());
    LinAlg::hessenberg(boltzmann2.data(), eedf2.data(), grid2.nCells());
    test_expr(eedf1.isApprox(eedf2));

    test_report;
    return nerrors;
}
