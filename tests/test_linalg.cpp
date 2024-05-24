/** \file
 *  Unit tests of the code in LinearAlgebra.h
 *
 *  \author Jan van Dijk
 *  \date   May 2024
 */

#include "LoKI-B/LinearAlgebra.h"
#include "tests/TestUtilities.h"
#include <algorithm>

using Matrix = loki::Matrix;

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

int main()
{
    testBandwidth();

    test_report;

    return nerrors;
}
