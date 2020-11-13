#ifndef LOKI_CPP_LINEARALGEBRA_H
#define LOKI_CPP_LINEARALGEBRA_H

//
// Created by daan on 15-5-19.
//

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

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace loki {

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

} // namespace loki

namespace loki {
namespace LinAlg {

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

double *hessenberg(const double *A, double *b, uint32_t n);
double *hessenberg(const double *A, double *b, uint32_t n, HessenbergWorkspace& hws);

} // namespace LinAlg
} // namespace loki

#endif // LOKI_CPP_LINEARALGEBRA_H

