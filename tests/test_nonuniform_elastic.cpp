/** \file
 *
 *  Unit tests for nonuniform elastic operator
 *
 *  \author Jan van Dijk
 *  \date   September 2023
 */

#include "LoKI-B/Grid.h"
#include "source/Operators.cpp"

#include "tests/TestUtilities.h"

int main()
{
    using namespace loki;

    Grid grid1(5,5);
    grid1.updateMaxEnergy(10);

    Grid::Vector ls(6); ls << 0.0,0.2,0.4,0.6,0.8,1.0;
    Grid grid2(ls,10,false); 

    Vector elasticCrossSection = Vector::Ones(6);
    ElasticOperator elasticOperator;
    SparseMatrix M1(5,5);
    SparseMatrix M2(5,5);
    elasticOperator.evaluate(grid1, elasticCrossSection, 300, M1);
    elasticOperator.evaluate(grid2, elasticCrossSection, 300, M2);
    test_expr(M1.isApprox(M2));

    test_report;
    return nerrors;
}
