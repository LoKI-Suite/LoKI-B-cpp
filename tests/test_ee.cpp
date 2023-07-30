/** \file
 *
 *  Unit tests for ElectronElectronOperator
 *
 *  \author Jan van Dijk
 *  \date   July 2023
 */

#include "LoKI-B/Grid.h"
#include "LoKI-B/Operators.h"
#include "LoKI-B/EedfUtilities.h"
#include "tests/TestUtilities.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>

double relativeError(double v1, double v2)
{
    return v1==0 && v2==0
        ? 0
        : std::abs(v1-v2)/(std::max(std::abs(v1),std::abs(v2)));
}

struct FluxTerms
{
    double up;
    double down;
    double net() const { return up - down; }
    double ref() const { return (std::abs(up)+std::abs(down))/2; }
    double relativeMismatch() const
    {
        return relativeError(up,down);
    } 
};

FluxTerms calculateFlux(const loki::ElectronElectronOperator& eeOp, const loki::Grid& grid, const loki::Vector& eedf, loki::Grid::Index n)
{
    /* For n==0 and n==grid.getNodes().size()-1, the flux is 0.
     * For other n, we calculate it from the coefficients A and B according to
     * G(node n) = N*gamma*du*( B_nf_n - A_(n-1)f_(n-1) ).
     */
    assert(n>=0 && n<grid.getNodes().size());
    if (n==0 || n==grid.getNodes().size()-1)
    {
        return {0,0};
    }
    else
    {
        const double aux = -loki::SI::gamma*grid.du();
        return { aux*eeOp.A[n-1]*eedf[n-1], aux*eeOp.B[n]*eedf[n] };
    }
}

/// \todo Also add non-Maxwellian tests as well, compare results with known refs.
void test1()
{
    using namespace loki;

    const unsigned nCells = 400;
    const double uMax = 4; // eV
    const double ne = 1e16; // m^-3
    const double n0 = 1e22; // m^-3

    const double shape = 1.0; // 1: Maxwell, 2: Druyvesteyn
    const double TeV = 1.0; // eV

    Grid grid(nCells,uMax);
    ElectronElectronOperator eeOperator;
    eeOperator.initialize(grid);
    eeOperator.updateABMatrices(grid);

    const Vector eedf = makePrescribedEDF(grid,shape,TeV);
    eeOperator.update_g_ee_AB(grid,eedf,ne,n0);

    std::vector<FluxTerms> flux(grid.getNodes().size());
    for (Grid::Index n = 0; n != grid.getNodes().size(); ++n)
    {
        flux[n] = calculateFlux(eeOperator,grid,eedf,n);
    }
    std::vector<double> src(grid.getCells().size());
    for (Grid::Index k = 0; k != grid.getCells().size(); ++k)
    {
        src[k] = -(flux[k+1].net()-flux[k].net())/grid.du();
    }

    std::cout << "#Node properties:" << std::endl;
    std::cout << "#energy\tup\tdown\tnet\trel.difference:" << std::endl;
    for (Grid::Index n = 0; n != grid.getNodes().size(); ++n)
    {
        std::cout << grid.getNodes()[n]
            << '\t' << flux[n].up
            << '\t' << flux[n].down
            << '\t' << flux[n].net()
            << '\t' << flux[n].relativeMismatch()
            << std::endl;
        test_expr( flux[n].relativeMismatch() < grid.getNodes().size()*std::numeric_limits<double>::epsilon() );
    }

    std::cout << std::endl << std::endl;
    std::cout << "#Cell properties: " << std::endl;
    std::cout << "#energy\tA/g\tB/g\tsrc:" << std::endl;
    double src_int = 0.0;
    for (Grid::Index k = 0; k != grid.nCells(); ++k)
    {
        /* for a Maxwellian, the correct value is 0.
         * Use the reference values that are provided by the fluxes to
         * obtain a meaningful reference to quantify the source error.
         */
        const double srcRef = (flux[k+1].ref()+flux[k].ref())/2/grid.du();
        const double srcRelError = std::abs(src[k])/srcRef;
        std::cout << grid.getCells()[k]
            << '\t' << eeOperator.A(k)/eeOperator.g_ee
            << '\t' << eeOperator.B(k)/eeOperator.g_ee
            << '\t' << src[k]
            << '\t' << srcRelError
            << std::endl;
        src_int += src[k]*grid.du();
        test_expr( srcRelError < grid.getCells().size()*std::numeric_limits<double>::epsilon() );
    }
    std::cout << "# src_int = " << src_int << std::endl;

    std::cout << std::endl << std::endl;
    double pElectronElectron;
    eeOperator.evaluatePower(grid,eedf,pElectronElectron);
    /** \todo We have power = (-SI::gamma * grid.du() * grid.du()) * (A - B).dot(eedf);
     * Calculate up and downfluxes separately, use those to obtain a reference value and
     * use that.
     */
    const double Eup = -SI::gamma * grid.du() * grid.du() * eeOperator.A.dot(eedf);
    const double Edn = -SI::gamma * grid.du() * grid.du() * eeOperator.B.dot(eedf);
    std::cout << "# up: " << Eup << ", down: " << Edn << ", net: " << (Eup-Edn) << std::endl;
    std::cout << "# pElectronElectron = " << pElectronElectron << std::endl;
}

int main()
{
    test1();

    test_report;
    return nerrors;
}
