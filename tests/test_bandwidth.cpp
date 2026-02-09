#include "LoKI-B/LinearAlgebra.h"
#include "tests/TestUtilities.h"

void test_diag()
{
    loki::Matrix M(4,4);
    M << 1.0,  0.0,  0.0, 0.0,
         0.0,  1.0,  0.0, 0.0,
         0.0,  0.0,  1.0, 0.0,
         0.0,  0.0,  0.0, 1.0;
    const loki::Bandwidth bw = loki::calculateBandwidth(M);
    test_expr(bw.first==0);
    test_expr(bw.second==0);
}

void test_tridiag()
{
    loki::Matrix M(4,4);
    M << 1.0,  1.0,  0.0, 0.0,
         1.0,  1.0,  1.0, 0.0,
         0.0,  1.0,  1.0, 1.0,
         0.0,  0.0,  1.0, 1.0;
    const loki::Bandwidth bw = loki::calculateBandwidth(M);
    test_expr(bw.first==-1);
    test_expr(bw.second==1);
}

void test_sltriangular()
{
    loki::Matrix M(4,5);
    M << 0.0,  0.0,  0.0, 0.0, 0.0,
         1.0,  0.0,  0.0, 0.0, 0.0,
         1.0,  1.0,  0.0, 0.0, 0.0,
         1.0,  1.0,  1.0, 0.0, 0.0;
    const loki::Bandwidth bw = loki::calculateBandwidth(M);
    test_expr(bw.first==-3);
    test_expr(bw.second==-1);
}

void test_uhessenberg()
{
    loki::Matrix M(5,5);
    M << 1.0,  1.0,  1.0, 1.0, 1.0,
         1.0,  1.0,  1.0, 1.0, 1.0,
         0.0,  1.0,  1.0, 1.0, 1.0,
         0.0,  0.0,  1.0, 1.0, 1.0,
         0.0,  0.0,  0.0, 1.0, 1.0;
    const loki::Bandwidth bw = loki::calculateBandwidth(M);
    test_expr(bw.first==-1);
    test_expr(bw.second==4);
}

/* The O matrix is a special case. The lower and upper bandwidths
 * are equal to +/- max(rows,cols).
 */
void test_null()
{
    loki::Matrix M(2,3);
    M << 0.0,  0.0,
         0.0,  0.0;
    const loki::Bandwidth bw = loki::calculateBandwidth(M);
    test_expr(bw.first==+3);
    test_expr(bw.second==-3);
}

int main()
{
    test_diag();
    test_tridiag();
    test_sltriangular();
    test_uhessenberg();
    test_null();

    test_report;

    return nerrors;
}
