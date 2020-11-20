/** \file
 *
 *  Interface of LoKI-B's linear algebra facilities.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2020 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \date   15. May 2019
 */

#ifndef LOKI_CPP_LINEARALGEBRA_H
#define LOKI_CPP_LINEARALGEBRA_H

/** \todo Should this define be done in CMakeLists?
 *  \todo Does USE_OPENMP have a meaning for the compiler? Or is it just
 *        used by LoKI-B to enable/disable OpenMP statements. In the latter
 *        case, it may be more clear to use LOKIB_USE_OPENMP.
 *  \todo If enabling OpenMP is controlled by the makefile (CMakeLists),
 *        USE_OPENMP is not needed, the #pragama statements can be made
 *        unconditional. Depending on the flags that are passed to the compiler,
 *        OpenMP will then be active (e.g. -fopenmp for gcc).
 */
#define USE_OPENMP

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

namespace loki {

/** In LoKI-B, the dense vector class of the Eigen library has been adopted
 *  as numerical vector type. For details on the properties of this class
 *  and its interfaces, we refer to the eigen documentation, see the
 *  documentation at https://eigen.tuxfamily.org
 */
using Vector = Eigen::VectorXd;

/** In LoKI-B, the dense matrix class of the Eigen library has been adopted
 *  as matrix type. See the documentation of Vector for a further discussion.
 */
using Matrix = Eigen::MatrixXd;

/** In LoKI-B, the SparseMatrix class of the Eigen library has been adopted
 *  as sparse matrix type. See the documentation of Vector for a further
 *  discussion.
 */
using SparseMatrix = Eigen::SparseMatrix<double>;

} // namespace loki

namespace loki {
namespace LinAlg {

/** The function hessenberg() requires temporary storage for vectors c, s and v
 *  of length n, the size of the solution vector. This storage is provided by
 *  an object of type HessenbergWorkspace.
 *
 *  \sa hessenberg
 *
 *  \author Jan van Dijk
 *  \date   November 2020
 */
struct HessenbergWorkspace
{
    HessenbergWorkspace(Vector::Index n)
     : c(n), s(n), v(n)
    {
    }
    Vector::Index size() const { return c.size(); }
    void resize(Vector::Index n)
    {
        c.resize(n);
        s.resize(n);
        v.resize(n);
    }
    Vector c;
    Vector s;
    Vector v;
};

/** Solve the matrix-vector equation Ax=b for an upper-Hessenberg matrix.
 *  Argument \a A is a pointer a dense matrix and must use column-major
 *  storage order. On input, \a b must point to the first element of the
 *  right-hand side. On return, \a b has been overwritten with the solution
 *  \a x of the problem. This pointer is also returned by the function.
 *  Argument \a n represents the lengths of \a b and x and the numbers of
 *  rows and columns of the (square) matrix \a A. It is up to the caller
 *  to ensure that \a A and \a b have these dimensions, and that \a A is
 *  indeed an upper-Hessenberg matrix.
 *
 *  \author Daan Boer
 *  \date   15. May 2019
 */
double *hessenberg(const double *A, double *b, uint32_t n);

/** An overload of hessenberg() that allows the caller to provide a workspace,
 *  preventing the allocation of temporary storage withing the hessenberg
 *  function. The size of workspace \a hws must be at least \a n, an error is
 *  generated if this is not the case. See the other overload of hessenberg
 *  for a discussion of the other arguments.
 *
 *  \sa HessenbergWorkspace
 *
 *  \author Jan van Dijk
 *  \date   November 2020
 */
double *hessenberg(const double *A, double *b, uint32_t n, HessenbergWorkspace& hws);

} // namespace LinAlg
} // namespace loki

#endif // LOKI_CPP_LINEARALGEBRA_H

