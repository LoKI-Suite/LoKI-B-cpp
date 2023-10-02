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

void write_eedf(std::ostream& os, const loki::Grid& grid, const loki::Vector& eedf, int uMax, double T, double EoN, double WoN)
{
    os << "# " << "nCells = " << grid.nCells() << std::endl;
    os << "# " << "uMax = " << uMax << std::endl;
    os << "# " << "T = " << T << std::endl;
    os << "# " << "EoN = " << EoN << std::endl;
    os << "# " << "WoN = " << WoN << std::endl;
    for (loki::Grid::Index k = 0; k < grid.nCells(); ++k)
    {
        os << grid.getCells()(k) << '\t' << eedf(k) << std::endl;
    }
}

int main()
{
    using namespace loki;

    const unsigned nCells = 1000;
    const double uMax = 2; // eV
    const double T = 0;
    const double eon = 1;
    const double won = 1;

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


    field1.row(0) = grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells());
    field2.row(0) = grid2.getCells().cwiseSqrt()*grid2.du();
    elastic1.row(0) = grid1.getCells().cwiseSqrt().cwiseProduct(grid1.duCells());
    elastic2.row(0) = grid2.getCells().cwiseSqrt()*grid2.du();

    Matrix total1 = field1 + elastic1;

    LinAlg::hessenberg(total1.data(), eedf1.data(), grid1.nCells());

    Vector eedfDruyvesteyn = Vector::Zero(nCells);
    double kT = Constant::kBeV*T;
    eedfDruyvesteyn << makePrescribedEDF(grid1,2,kT);
 
    LinAlg::hessenberg(field1.data(), eedf1.data(), grid1.nCells());
    // writeGnuplot(std::cout, "Eedf1", "Energy (eV)", "Eedf (eV^-3/2)", grid1.getCells(), eedf1);

    // std::cout << "set multiplot" << std::endl;
    // writeGnuplot(std::cout, "Eedf1", "Energy (eV)", "Eedf (eV^-3/2)", grid1.getCells(), eedf1);
    // writeGnuplot(std::cout, "Eedf2", "Energy (eV)", "Eedf (eV^-3/2)", grid1.getCells(), eedfDruyvesteyn);
    // std::cout << "unset multiplot" << std::endl;

    std::ofstream ofs("eedfDruyvesteyn.dat");
    write_eedf(ofs, grid1, eedf1, uMax, T, eon, won);

    test_report;
    return nerrors;
}
