/** \file
 *
 *  Unit tests for nonuniform CAR operator
 *
 *  \author Jop Hendrikx
 *  \date   October 2023
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/EedfUtilities.h"
#include "tests/TestUtilities.h"
#include "source/Operators.cpp"

int main()
{
    using namespace loki;
    
    int nCells = 1000;
    double uMax = 2;
    double Tg = 300;
    double kT = Constant::kBeV*Tg;
    
    auto gas = Gas("carGas");
    gas.anharmonicFrequency = 1;
    gas.electricQuadrupoleMoment = 1;
    gas.harmonicFrequency = 1;
    gas.rotationalConstant = 1;
    gas.fraction = 1;

    std::vector<const loki::Gas*> carGases;
    carGases.push_back(&gas);
 
    Grid::Vector nodeDistribution = Vector::LinSpaced(nCells + 1, 0.0, 1.0);
    Grid grid1(nodeDistribution, uMax, false);
    Grid grid2(nCells, uMax);

    CAROperator caroperator(carGases);

    SparseMatrix M1(nCells, nCells);
    SparseMatrix M2(nCells, nCells);

    caroperator.evaluate(grid1, Tg, M1);
    caroperator.evaluate(grid2, Tg, M2);
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

    Vector eedfMaxwell = Vector::Zero(nCells);
    eedfMaxwell << makePrescribedEDF(grid1,1,kT);

    // Calculate relative error
    Vector relativeError = (eedf1 - eedfMaxwell).cwiseQuotient(eedfMaxwell).cwiseAbs();
    test_expr((relativeError.mean() <  0.02));

    test_report;
    return nerrors;
}
