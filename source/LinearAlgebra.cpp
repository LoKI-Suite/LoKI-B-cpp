/** \file
 *
 *  Implementation of LoKI-B's linear algebra facilities.
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
 *  \date   15 May 2019
 */

#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Log.h"

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <type_traits>

namespace loki {

double maxRelDiff(const Vector& v1, const Vector& v2)
{
//#define LOKIB_USE_OLD_MAXRELDIFF
#ifdef LOKIB_USE_OLD_MAXRELDIFF
    const double res = ((v2 - v1).cwiseAbs().array() / v2.array()).maxCoeff();
#else
    double res=-1.0;
    assert(v1.size()==v2.size());
    for (Vector::Index i=0; i!=v1.size(); ++i)
    {
        const double ref = (std::abs(v1(i)) + std::abs(v2(i)))/2;
        if (ref!=0)
        {
            const double dif = std::abs(v1(i)-v2(i));
            res = std::max(res,dif/ref);
        }
    }
    if (res<0.0)
    {
        throw std::runtime_error("Maximum relative difference: vectors empty or zero.");
    }
#endif
    // std::cout << "maxRelDiff: " << res << std::endl;
    return res;
}

std::pair<Matrix::Index,Matrix::Index> calculateBandwidth(const Matrix& m)
{
    /** \todo In principle, this can be optimized, since for each row the
     *  columns that are within bands that are in [min,max] do not need to
     *  be checked another time.
     */
    static_assert(std::is_signed_v<Matrix::Index>);
    /* 'illegal' starting values. For a zero-valued matrix
     * these will also be the final values.
     */
    Matrix::Index min=+std::max(m.rows(),m.cols());
    Matrix::Index max=-std::max(m.rows(),m.cols());
    for (Matrix::Index r=0; r!=m.rows(); ++r)
    for (Matrix::Index c=0; c!=m.cols(); ++c)
    {
        if (m(r,c)!=0.0)
        {
            const Matrix::Index del = c-r;
            min = std::min(min,del);
            max = std::max(max,del);
        }
    }
    return { min, max };
}

} // namespace loki

namespace loki {
namespace LinAlg {

namespace impl {

void givens(double a, double b, double &c, double &s)
{
    const double root = std::hypot(a,b);

    c =  a / root;
    s = -b / root;
}

} // namespace impl

double *hessenberg(const double *A, double *b, uint32_t n, HessenbergWorkspace& hws)
{
    if (hws.size()<n)
    {
        Log<Message>::Error("hessenberg: workspace is too small.");
    }
    auto *c = hws.c.data();
    auto *s = hws.s.data();
    auto *v = hws.v.data();
    double *x = b;

    // Matrix A is column major
    std::memcpy(v, A + (n - 1) * n, n * sizeof(double));

    for (uint32_t k = n - 1; k > 0; --k)
    {
        impl::givens(v[k], A[(k - 1) * n + k], c[k], s[k]);

        x[k] /= c[k] * v[k] - s[k] * A[(k - 1) * n + k];

        for (int64_t i = 0; i < k - 1; ++i)
        {
            x[i] += x[k] * s[k] * A[(k - 1) * n + i] - x[k] * c[k] * v[i];
            v[i] = c[k] * A[(k - 1) * n + i] + s[k] * v[i];
        }

        x[k - 1] += x[k] * s[k] * A[(k - 1) * n + k - 1] - x[k] * c[k] * v[k - 1];
        v[k - 1] = c[k] * A[(k - 1) * n + k - 1] + s[k] * v[k - 1];
    }

    double t1 = x[0] / v[0];

    for (uint32_t k = 1; k < n; ++k)
    {
        const double t2 = x[k];
        x[k - 1] = c[k] * t1 - s[k] * t2;
        t1 = c[k] * t2 + s[k] * t1;
    }

    x[n - 1] = t1;

    return x;
}

double *hessenberg(const double *A, double *b, uint32_t n)
{
    HessenbergWorkspace hws(n);
    return hessenberg(A,b,n,hws);
}

void solveTDMA(const Matrix& A, const Vector& b, Vector& x)
{
    static_assert(std::is_signed_v<Vector::Index>);
    Vector P( x.size() );
    Vector Q( x.size() );

    const Vector::Index last=x.size()-1;

    // 1. forward sweep. Calculares auxiliary vectors P and Q from A and b.

    if (A(0,0)==0.0)
    {
            throw std::runtime_error("Solver error: ap==0.");
    }
    P[0] = A(0,1)/A(0,0);
    Q[0] =   b[0]/A(0,0);
    for (Vector::Index p=1; p <= last ; ++p)
    {
            const typename Matrix::Scalar denom = A(p,p)-A(p,p-1)*P[p-1];
            if (denom==0.0)
            {
                     throw std::runtime_error("TDMA solver error: zero denominator.");
            }
            const typename Matrix::Scalar _1_denom = 1.0/denom;
            P[p] = (p==last) ? 0.0 : A(p,p+1)*_1_denom;
            Q[p] = (b[p]-A(p,p-1)*Q[p-1])*_1_denom;
    }

    // 2. backward sweep. Calculates x from P and Q.

    x[last]=Q[last];
    for (Vector::Index p=last-1; p>=0; --p)
    {
            x[p] = Q[p] - P[p]*x[p+1];
    }
}

void solveTDMA(const Matrix& A, Vector& bx)
{
    solveTDMA(A,bx,bx);
}

} // namespace LinAlg
} // namespace loki
