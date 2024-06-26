/** \file
 *  Unit tests of the code in LinearAlgebra.h
 *
 *  \author Jan van Dijk
 *  \date   May 2024
 */

#include "LoKI-B/LinearAlgebra.h"
#include "tests/TestUtilities.h"
#include <algorithm>
#include <chrono>
#include <limits>

using Matrix = loki::Matrix;
using Vector = loki::Vector;

void doTestBandwidth(Matrix::Index exp_min, Matrix::Index exp_max, const Matrix& m)
{
    const auto bw = loki::calculateBandwidth(m);
    test_expr(bw.first == exp_min);
    test_expr(bw.second == exp_max);
    const auto bwT = loki::calculateBandwidth(m.transpose());
    test_expr(bwT.first == -exp_max);
    test_expr(bwT.second == -exp_min);
}

void testBandwidth()
{
    doTestBandwidth(  2, -2, (Matrix(2,2) << 0.0, 0.0, 0.0, 0.0 ).finished() );
    doTestBandwidth(  0,  0, (Matrix(2,2) << 1.0, 0.0, 0.0, 1.0 ).finished() );
    doTestBandwidth(  0,  1, (Matrix(2,2) << 1.0, 1.0, 0.0, 1.0 ).finished() );
    doTestBandwidth(  0,  1, (Matrix(3,2) << 1.0, 1.0, 0.0, 1.0, 0.0, 0.0).finished() );
    doTestBandwidth( -1, -1, (Matrix(2,3) << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0).finished() );
}

void solveLU(const Matrix& A, const Vector& b, Vector& x)
{
    x = A.partialPivLu().solve(b);
}

void solveHessenberg(const Matrix& A, const Vector& b, Vector& x)
{
    x = b;
    loki::LinAlg::hessenberg(A.data(),x.data(),x.size());
}

void solveTDMA(const Matrix& A, const Vector& b, Vector& x)
{
    loki::LinAlg::solveTDMA(A,b,x);
}

void solveMatrix(const Matrix& A, const Vector& b, Vector& x, const Vector& x_anal, const std::string& name, void(*solver)(const Matrix&, const Vector&, Vector&))
{
    using namespace std::chrono;
    auto begin = high_resolution_clock::now();
    (*solver)(A,b,x);
    auto end = high_resolution_clock::now();
    std::cout << "Name: " << name << ", duration: " << duration_cast<microseconds>(end - begin).count() << "mus" << std::endl;
    test_expr((x-x_anal).cwiseQuotient(x_anal).cwiseAbs().sum()/x.size()<100*x.size()*std::numeric_limits<double>::epsilon());
}

/** Set up, test and benchmark solving a tridiagonal system (a diffusion
 *  equation) with the available solvers.
 */
void testSolversTridiagonal(Matrix::Index n)
{
    /* set up a diffusion equation with diffusion coefficient 1, subject to
     * boundary conditions x(0)=1, x(1)=2. The analytical solution is given
     * by x_anal(c)=1+c.
     */
    const double dx = 1.0/(n-1);
    Vector x_anal(n);
    for (Matrix::Index i=0; i!=n; ++i)
    {
        x_anal[i] = 1.0+i*dx;
    }
    Matrix A(n,n);
    A.setZero();
    A.diagonal(-1).fill(-1.0/dx/dx);
    A.diagonal( 0).fill(+2.0/dx/dx);
    A.diagonal(+1).fill(-1.0/dx/dx);
    Vector b(n);
    b.setZero();
    A(0,0)=1.0/dx/dx;
    A(0,1)=0.0;
    b(0)=1/dx/dx;
    A(n-1,n-1)=1.0/dx/dx;
    A(n-1,n-2)=0.0;
    b(n-1)=2.0/dx/dx;

    Vector x(n);
    x.setZero();
    solveMatrix(A,b,x,x_anal,"LU",&solveLU);
    solveMatrix(A,b,x,x_anal,"Hessenberg",&solveHessenberg);
    solveMatrix(A,b,x,x_anal,"TDMA",&solveTDMA);
}

int main()
{
    testBandwidth();
    testSolversTridiagonal(1001);

    test_report;

    return nerrors;
}
