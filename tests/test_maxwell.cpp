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
    const double uMax = 10; // eV

    Vector elasticCrossSection = Vector::Ones(nCells+1);

    Grid::Vector ls(nCells+1);
    ls << Vector::LinSpaced(nCells + 1, 0.0, 1.0);
    Grid grid1(ls,uMax,false);
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
    Matrix boltzmann2 = M2.toDense();
    boltzmann1.row(0) = grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells());
    boltzmann2.row(0) = grid2.getCells().cwiseSqrt()*grid2.du();
    test_expr(boltzmann1.isApprox(boltzmann2));

    LinAlg::hessenberg(boltzmann1.data(), eedf1.data(), grid1.nCells());
    LinAlg::hessenberg(boltzmann2.data(), eedf2.data(), grid2.nCells());
    test_expr(eedf1.isApprox(eedf2));

    // Analytical solution
    double kT = Constant::kBeV*T;
    Vector eedfMaxwell = makePrescribedEDF(grid1,1,kT);

    // Calculate relative error
    Vector relativeError = (eedf1 - eedfMaxwell).cwiseAbs().cwiseQuotient(eedfMaxwell.cwiseAbs());
    std::cerr << relativeError.mean() << std::endl;
    test_expr(eedf1.isApprox(eedfMaxwell));

    // Plot nonuniform numerical and analytical solution
    // std::cout << "set multiplot" << std::endl;
    // writeGnuplot(std::cout, "Eedf1", "Energy (eV)", "Eedf (eV^-3/2)", grid1.getCells(), eedf1);
    // writeGnuplot(std::cout, "Eedf2", "Energy (eV)", "Eedf (eV^-3/2)", grid1.getCells(), eedfMaxwell);
    // std::cout << "unset multiplot" << std::endl;

    // Plot relative error
    // writeGnuplot(std::cout, "Eedf1", "Energy (eV)", "MRE", grid1.getCells(), relativeError);

    test_report;
    return nerrors;
}
