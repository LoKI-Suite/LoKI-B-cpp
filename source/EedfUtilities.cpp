/** \file
 *
 *  Implementation of LoKI-B EEDF-related utility functions.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk (C++ version)
 *  \date   29 July 2023
 */

#include "LoKI-B/EedfUtilities.h"
#include "LoKI-B/GridOps.h"
#include <cmath>

namespace loki {

double getMeanEnergy(const Vector& edf, const Grid& grid)
{
    return energyIntegral(grid,grid.getCells().cwiseProduct(grid.getCells().cwiseSqrt()),edf);
}

void normalizeEDF(Vector& edf, const Grid& grid)
{
    edf /= energyIntegral(grid,grid.getCells().cwiseSqrt(),edf);
}

void makePrescribedEDF(Vector& edf, const Grid& grid, double s, double T_eV, bool normalize)
{
    edf.resize(grid.nCells());
    const double gamma_3_2s = std::tgamma(3/(2*s));
    const double gamma_5_2s = std::tgamma(5/(2*s));
    const double factor = s*std::pow(gamma_5_2s,1.5) * std::pow(gamma_3_2s,-2.5) * std::pow(2/(3*T_eV),1.5);
    for (Grid::Index k = 0; k != grid.nCells(); ++k)
    {
        const double energy = grid.getCell(k);
        edf(k) = factor * std::exp(-std::pow(energy*gamma_5_2s*2/(gamma_3_2s*3*T_eV),s));
    }
    if (normalize)
    {
        normalizeEDF(edf,grid);
    }
}

Vector makePrescribedEDF(const Grid& grid, double s, double T_eV, bool normalize)
{
    Vector edf(grid.nCells());
    makePrescribedEDF(edf,grid,s,T_eV,normalize);
    return edf;
}

void solveEEDF(Vector &eedf, Matrix &matrix, const Grid &grid)
{
// show matrix inversion times?
#define LOKIB_MATRIX_TIME_INVERSION 0

// print the bandwidth and solver choice to the console?
#define LOKIB_MATRIX_SHOW_BANDWIDTHS 0

#if LOKIB_MATRIX_TIME_INVERSION
    using namespace std::chrono;
    const auto begin = high_resolution_clock::now();
#endif

    /* 1. In this function we use the matrix, the grid and the eedf vector.
     *    First check that the dimensions of these variables ae consistent.
     */
    if (matrix.rows()!=matrix.cols())
    {
        throw std::runtime_error("invertMatrix: the matrix is not square.");
    }
    if (matrix.rows()!=grid.nCells())
    {
        throw std::runtime_error("invertMatrix: the matrix dimensions do not match the grid's number of cells.");
    }
    if (matrix.rows()!=eedf.size())
    {
        throw std::runtime_error("invertMatrix: the matrix dimensions do not match the length of the solution vector.");
    }
    /* In principle we can allow a matrix with a single row and column, but that
     * is not of practical interest, and in the code below we use the element
     * matrix(1,1).
     */
    if (matrix.rows()<2)
    {
        throw std::runtime_error("invertMatrix: the matrix must have at least 2 rows.");
    }

    /* 2. Make the system non-singular by fixing the first value of the EEDF to 1:
     *   replace the first equation with M(1,1)*eedf[0] = M(1,1). After solving
     *   the system, the eedf is rescaled as to satisfy the normalization
     *   condition. Choosing M(1,1) is semi-arbitrary, but avoids that we
     *   unnecessarily increase the dynamic range of the matrix elements. In
     *   principle any other non-zero value will do.
     */
    matrix.row(0).setZero();
    matrix(0,0)=matrix(1,1);
    eedf.setZero();
    eedf[0] = matrix(1,1);

    /* 3. Calculate the lower and upper bandwidth and use the most efficient
     *    algorithm for solving the matrix. For a tridiagonal matrix (with lower
     *    and upper band widths [-1,1]) use TDMA, for an upper Hessenberg matrix
     *    (with bandwidths [-1,X]), use the hessenberg hessenberg solver.
     *    Otherwise, use eigen's LU solver with partial pivoting.
     */
#if LOKIB_MATRIX_TIME_INVERSION
    const auto begin_bw = high_resolution_clock::now();
#endif
    const auto bw = calculateBandwidth(matrix);
#if LOKIB_MATRIX_TIME_INVERSION
    const auto end_bw = high_resolution_clock::now();
#endif
#if LOKIB_MATRIX_SHOW_BANDWIDTHS==1
    std::cout << "In ElectronKineticsBoltzmann::invertMatrix" << std::endl;
    std::cout << " * bandwidth: [" << bw.first << ',' << bw.second << ']' << std::endl;
#endif
    if (bw.first==-1 && bw.second==+1)
    {
        LinAlg::solveTDMA(matrix,eedf);
    }
    else if (bw.first==-1)
    {
        /* solve Ax=b. On entry of LinAlg::hessenberg, the second argument of
         * LinAlg::hessenberg (eedf.data) must point to b. After returning,
         * this vector has been overwritten with the solution vector x of the
         * system of equations.
         */
        LinAlg::hessenberg(matrix.data(), eedf.data(), grid.nCells());
    }
    else
    {
        // at this point, eedf has been set up to contain b.
        Vector b(eedf);
        eedf = matrix.partialPivLu().solve(b);
    }
    /* 4. We have calculated the new EEDF, but eedf[0]==1 and the result is
     *    correct up to a multiplicative constant. Scale the EEDF to satisfy
     *    the normalization condition (the integral of f(u)sqrt(u) over the
     *    energies must be unity) to make it final.
     */
    normalizeEDF(eedf,grid);

#if LOKIB_MATRIX_TIME_INVERSION
    const auto end = high_resolution_clock::now();
    std::cerr << " invertMatrix time report:"
        << "\n * bandwidth: " << duration_cast<microseconds>(end_bw - begin_bw).count() << "mus"
        << "\n * total:     " << duration_cast<microseconds>(end - begin).count() << "mus"
        << std::endl;
#endif
}



} // namespace loki

