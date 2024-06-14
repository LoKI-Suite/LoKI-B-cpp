/** \file
 *
 *  Declarations of functions that implement grid operations (interpolation,
 *  differentiation, integration).
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

#ifndef LOKI_CPP_GRIDOPS_H
#define LOKI_CPP_GRIDOPS_H

#include "LoKI-B/Grid.h"
#include <stdexcept>

namespace loki {

    /** Interpolate the values \a src that are defined on the faces of \a grid to the
     *  cells, store the result in \a tgt. Argument \a tgt is resized if necessary.
     *
     *  If the size of \a src does not match the number of faces of \a grid, a
     *  std::runtime_error is thrown.
     *
     *  \author Jan van Dijk
     *  \date   13 June 2024
     */
    inline void interpolateNodalToCell(const Grid& grid, const Vector& src, Vector& tgt)
    {
        if (src.size()!=grid.getNodes().size())
        {
            throw std::runtime_error("interpolateNodalToCell: grid and source vector have incompatible sizes.");
        }
        tgt = (src.head(grid.nCells()) + src.tail(grid.nCells()))/2;
    }

    /** Interpolate the values \a src that are defined on the faces of \a grid to the
     *  cells, return the result. This is implemented by calling the three-argument
     *  overload of this function on a local vector, and returning the result. See that
     *  function for more information.
     *
     *  \author Jan van Dijk
     *  \date   13 June 2024
     */
    inline Vector interpolateNodalToCell(const Grid& grid, const Vector& src)
    {
        Vector tgt;
        interpolateNodalToCell(grid,src,tgt);
        return tgt;
    }

    /** Calculate and return an approximation of the derivative of the field
     *  \a f. Both the field \a f and the resulting vactor are defined in the
     *  cells of the \a grid.
     *
     *  \author Jan van Dijk
     *  \date   14 June 2024
     */
    inline void cellDerivative(Vector& fprime, const Grid& grid, const Vector& f)
    {
        const Grid::Index n = grid.nCells();
        if (f.size()!=n)
        {
            throw std::runtime_error("cellToNodalDerivative: "
                                     "field has wrong size.");
        }
        fprime.resize(n);
        if (grid.isUniform())
        {
            fprime[0] = (f[1] - f[0]) / grid.du();
            fprime[n - 1] = (f[n - 1] - f[n - 2]) / grid.du();
            fprime.segment(1, n - 2) = (f.segment(2, n - 2) - f.segment(0, n - 2)) / (2 * grid.du());
        }
        else
        {
            fprime[0] = (f[1] - f[0]) / grid.duNode(1);
            fprime[n - 1] = (f[n - 1] - f[n - 2]) / grid.duNode(n-1);
            fprime.segment(1, n - 2) = (f.segment(2, n - 2) - f.segment(0, n - 2))
                .cwiseQuotient(grid.duNodes().segment(1, n - 2) + grid.duNodes().segment(2, n - 2));
        }
    }

    inline Vector cellDerivative(const Grid& grid, const Vector& f)
    {
        Vector fprime;
        cellDerivative(fprime,grid,f);
        return fprime;
    }

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
        if (N<2)
        {
            throw std::runtime_error("fgPrimeEnergyIntegral: "
                                     "grid must contain at least two cells.");
        }
        if (f.size()!=g.size() || f.size()!=N)
        {
            throw std::runtime_error("fgPrimeEnergyIntegral: "
                                     "multiplicands have wrong sizes.");
        }
        if (grid.isUniform())
        {
            /** \todo The parenthesized expressions below are the coefficients
             *  c such that dot(c,f) yields the result. Since we only care about
             *  the results, not about the coefficients as such, we may decide
             *  to also implement the isUniform() path in the simpler way that
             *  is also used for the non-uniform case.
             */
            // Note that the result is grid-independent in this case (as long as
            // it is uniform).
            double res = 0.0;
            switch (N)
            {
                // NOTE: N>=2 has been checked above
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
            // NOTE: N>=2 has been checked above
            double res = 0.0;
            res += f[  0]*(g[  1]-g[  0])/grid.duNode(1)*grid.duCell(0);
            res += f[N-1]*(g[N-1]-g[N-2])/grid.duNode(N-1)*grid.duCell(N-1);
            for (Grid::Index k=1; k<=N-2; ++k)
            {
                const double dum = grid.duNode(k);
                const double dup = grid.duNode(k+1);
                res += f[k]*(g[k+1] - g[k-1])/(dum+dup)*grid.duCell(k);
            }
            return res;
        }
    }

} // namespace loki

#endif // LOKI_CPP_GRIDOPS_H

