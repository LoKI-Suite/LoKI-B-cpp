/** \file
 *
 *  Unit tests for nonuniform field operator
 *
 *  \author Jop Hendrikx
 *  \date   September 2023
 */

#include "LoKI-B/Grid.h"
#include "source/Operators.cpp"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Gnuplot.h"
#include "LoKI-B/EedfUtilities.h"

#include "tests/TestUtilities.h"


int main()
{
    using namespace loki;


    const unsigned nCells = 1024;
    const double uMax = 2; // eV
    const double eon = 5;
    const double won = 0;

    Vector fieldCrossSection = Vector::Ones(nCells+1);
    Vector elasticCrossSection = Vector::Ones(nCells+1);

    auto nodeDistribution = Grid::Vector::LinSpaced(nCells + 1, 0.0, 1.0);
    Grid grid1(nodeDistribution,uMax,false);
    Grid grid2(nCells, uMax);

    FieldOperator fieldOperator1(grid1);
    FieldOperator fieldOperator2(grid2);
    SparseMatrix M1(nCells,nCells);
    SparseMatrix M2(nCells,nCells);

    const double CIEff = 0.0;
    fieldOperator1.evaluate(grid1, fieldCrossSection, eon, won, CIEff, M1);
    fieldOperator2.evaluate(grid2, fieldCrossSection, eon, won, CIEff, M2);
    test_expr(M1.isApprox(M2));

    Vector eedf1 = Vector::Zero(nCells);
    eedf1[0] = 1.;

    Vector eedf2 = Vector::Zero(nCells);
    eedf2[0] = 1.;

    Matrix field1 = M1.toDense();
    Matrix field2 = M2.toDense();

    field1.row(0) = grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells());
    field2.row(0) = grid2.getCells().cwiseSqrt()*grid2.du();

    LinAlg::hessenberg(field1.data(), eedf1.data(), grid1.nCells());
    LinAlg::hessenberg(field2.data(), eedf2.data(), grid2.nCells());

    test_expr( eedf1.isApprox(eedf2));
    double analytical = 3./2. / uMax / std::sqrt(uMax);
    /** \todo Even when eedf is constant, as it will be for the present case,
     *  it will be off by a constant factor because of the discretization error
     *  in the normalization constant. The computed normalized constant eedf
     *  will be equal to the commented out expression below.
     */
    // double analytical = 1./grid1.getCells().cwiseSqrt().dot(grid1.duCells());
    test_expr( (eedf1-Vector::Ones(nCells)*analytical).cwiseAbs().maxCoeff() / analytical < 0.003 );

    test_report;
    return nerrors;
}
