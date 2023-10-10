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

#include <iostream>
#include <fstream>

#include "tests/TestUtilities.h"


int main()
{
    using namespace loki;

    const unsigned nCells = 1000;
    const double uMax = 20; // eV
    const double T = 0;
    const double eon = 10;
    const double won = 0;

    Vector fieldCrossSection = Vector::Ones(nCells+1);
    Vector elasticCrossSection = Vector::Ones(nCells+1);

    Grid::Vector ls(nCells+1); 
    ls << Vector::LinSpaced(nCells + 1, 0.0, 1.0);
    Grid grid1(ls,uMax,false); 
    Grid grid2(nCells, uMax);

    FieldOperator fieldOperator1(grid1);
    FieldOperator fieldOperator2(grid2);
    SparseMatrix M1(nCells,nCells);
    SparseMatrix M2(nCells,nCells);
    
    fieldOperator1.evaluate(grid1, fieldCrossSection, eon, won, M1);
    fieldOperator2.evaluate(grid2, fieldCrossSection, eon, won, M2);

    ElasticOperator elasticOperator;
    SparseMatrix Melastic1(nCells, nCells);
    SparseMatrix Melastic2(nCells, nCells);
    elasticOperator.evaluate(grid1, elasticCrossSection, T, Melastic1);
    elasticOperator.evaluate(grid2, elasticCrossSection, T, Melastic2);


    Vector eedf1 = Vector::Zero(nCells);
    eedf1[0] = 1.;

    Vector eedf2 = Vector::Zero(nCells);
    eedf2[0] = 1.; 

    Matrix field1 = M1.toDense();
    Matrix field2 = M2.toDense();
    Matrix elastic1 = Melastic1.toDense();
    Matrix elastic2 = Melastic2.toDense();

    Matrix total1 = field1 + elastic1;
    total1.row(0) = grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells());

    LinAlg::hessenberg(total1.data(), eedf1.data(), grid1.nCells());
    Vector eedf = eedf1/(eedf1.dot(grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells())));

   

    double averageEnergy = eedf.dot(grid1.getCells().cwiseProduct(grid1.getCells().cwiseSqrt()).cwiseProduct(grid1.duCells()));
    Vector eedfDruyvesteyn = Vector::Zero(nCells);
    double kT = 2./3.*averageEnergy;

    eedfDruyvesteyn << makePrescribedEDF(grid1,2,kT);

    double relativeError = ((eedf1-eedfDruyvesteyn).cwiseQuotient(eedfDruyvesteyn).cwiseAbs()).mean();
    test_expr((relativeError <  0.005));

    test_report;
    return nerrors;
}
