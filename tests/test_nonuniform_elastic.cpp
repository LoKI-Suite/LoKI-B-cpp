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

    Grid::Vector ls(nCells+1); 
    ls << Vector::LinSpaced(nCells + 1, 0.0, 1.0);
    Grid grid1(ls,uMax,false); 
    Grid grid2(nCells, uMax);

    ElasticOperator elasticOperator;
    SparseMatrix M1(nCells,nCells);
    SparseMatrix M2(nCells,nCells);
    
    elasticOperator.evaluate(grid1, elasticCrossSection, T, M1);
    elasticOperator.evaluate(grid2, elasticCrossSection, T, M2);
    test_expr(M1.isApprox(M2));

    test_report;
    return nerrors;
}
