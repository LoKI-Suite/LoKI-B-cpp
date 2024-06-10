/** \file
 *
 *  Declarations of functions for evaluating integrals.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \author Jan van Dijk and Daan Boer
 *  \date   June 2024
 */

#ifndef LOKI_CPP_INTEGRALS_H
#define LOKI_CPP_INTEGRALS_H

#include "LoKI-B/Grid.h"
#include <stdexcept>

namespace loki {


    /** Calculate and return an approximation of \f$ \gamma\int_0^{u_max} f(u)g(u)du \f$.
     *  Arguments \a f and \a g must be defined in the cells of the \a grid.
     *  The function returns \f$ \sum_k f_kg_k(\Delta u)_k \f$; for uniform grids
     *  Delta u is constant and this reduces to \f$ (\sum_k f_kg_k)\Delta u \f$.
     *
     *  \author Jan van Dijk
     *  \date   5 June 2024
     */
    template <class VectorExpr>
    double energyIntegral(const Grid& grid, const VectorExpr& f, const Vector& g)
    {
        const Grid::Index N = grid.nCells();
        if (f.size()!=g.size() || f.size()!=N)
        {
            throw std::runtime_error("energyIntegral: "
                                     "multiplicands have wrong sizes.");
        }
        if (grid.isUniform())
        {
            return f.dot(g)*grid.du();
        }
        else
        {
            return f.cwiseProduct(g).dot(grid.duCells());
            throw std::runtime_error("energyIntegral: non-equidistant "
                                     "meshes are not yet supported.");
        }
    }

    /** Calculate coefficients \a C such that C.dot(f) is
     *  an approximation of \f$ \gamma\int_0^{u_max} g\frac{df}{du}du \f$.
     *  Argument \a C is resized to D0.size() by this function.
     *
     *  \todo Document assumptions on the grid layout.
     *
     *  \author Jan van Dijk
     *  \date   5 June 2024
     */
    template <class VectorExpr>
    void fgPrimeEnergyIntegralCoefficients(Vector& C, const VectorExpr& D0)
    {
        const Grid::Index N = D0.size();
        if (N<2)
        {
            throw std::runtime_error("doCalculateD0dfduIntCoefD0dfduIntCoefs: "
                                     "grid must contain at least two cells.");
        }
        C.resize(N);
        switch (N)
        {
            case 2:
                C[0] = -D0[0] - D0[1];
                C[1] = -C[0];
            break;
            case 3:
                C[0] = -D0[0]   - D0[1]/2;
                C[1] = +D0[0]   - D0[2];
                C[2] = +D0[1]/2 - D0[2];
            break;
                default:
                C[0] = -D0[0] - D0[1]/2;
                C[1] = +D0[0] - D0[2]/2;
                for (Grid::Index k=2; k<=N-3; ++k)
                {
                    C[k] = (D0[k-1] - D0[k+1])/2;
                }
                C[N-2] = +D0[N-3]/2 - D0[N-1];
                C[N-1] = +D0[N-2]/2 + D0[N-1];
            break;
        }
    }

    /** Calculate and return an approximation of \f$ \gamma\int_0^{u_max} f\frac{dg}{du}du \f$.
     *  Arguments \a f and \a g must be defined in the cells of the \a grid.
     *
     *  \author Jan van Dijk
     *  \date   5 June 2024
     */
    template <class VectorExpr>
    double fgPrimeEnergyIntegral(const Grid& grid, const VectorExpr& f, const Vector& g)
    {
        const Grid::Index N = grid.nCells();
        if (f.size()!=g.size() || f.size()!=N)
        {
            throw std::runtime_error("fgPrimeEnergyIntegral: "
                                     "multiplicands have wrong sizes.");
        }
        if (grid.isUniform())
        {
            // Note that the result is grid-independent in this case (as long as
            // it is uniform).
            double res = 0.0;
            switch (N)
            {
                case 0:
                case 1:
                    throw std::runtime_error("fgPrimeEnergyIntegral: "
                                             "grid must contain at least two cells.");
                break;
                case 2:
                    res += (-f[0] - f[1])*g[0];
                    res += ( f[0] + f[1])*g[1];
                break;
                case 3:
                    res += (-f[0]   - f[1]/2)*g[0];
                    res += (+f[0]   - f[2]  )*g[1];
                    res += (+f[1]/2 - f[2]  )*g[2];
                break;
                default:
                    res += (-f[0] - f[1]/2)*g[0];
                    res += (+f[0] - f[2]/2)*g[1];
                    for (Grid::Index k=2; k<=N-3; ++k)
                    {
                        res += (f[k-1] - f[k+1])/2*g[k];
                    }
                    res += (+f[N-3]/2 - f[N-1])*g[N-2];
                    res += (+f[N-2]/2 + f[N-1])*g[N-1];
                break;
            }
            return res;
        }
        else
        {
            throw std::runtime_error("fgPrimeEnergyIntegral: non-equidistant "
                                     "meshes are not yet supported.");
        }
    }

} // namespace loki

#endif // LOKI_CPP_INTEGRALS_H

