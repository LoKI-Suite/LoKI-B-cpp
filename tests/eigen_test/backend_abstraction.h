
#define EIGEN

#ifdef EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif

/* FUNCTIONS TO IMPLEMENT
 *
 * 1. Binary Expressions: cwiseProduct, cwiseQuotient, dot
 * 2. Unary Expressions: cwiseSqrt, cwiseAbs, cwiseAbs2
 * 
 * For Sparse Matrices
 * 
 * 1. Getters: coeffRef, coeff
 * 2. A way to make and set patterns
 * 
 */

namespace loki
{
namespace la
{
template <typename Library, typename Type>
struct Traits;

template <typename Library, typename type>
struct Functions;

#ifdef EIGEN
class L_Eigen;

template <typename Type>
struct Traits<L_Eigen, Type>
{
    using SparseMatrix = Eigen::SparseMatrix<Type>;
    using Matrix = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Type, Eigen::Dynamic, 1>;

    template <typename Derived>
    using AnyDense = Eigen::PlainObjectBase<Derived>;

    using LUSolveType = const Eigen::Solve<Eigen::PartialPivLU<Matrix>, Vector>;
};

template <typename Type>
struct Functions<L_Eigen, Type>
{
    using T = Traits<L_Eigen, Type>;

    static typename T::LUSolveType pLUSolve(const typename T::Matrix &A, const typename T::Vector &b)
    {
        return A.partialPivLu().solve(b);
    }

    template <typename Derived, typename... Args>
    static void setRandom(typename T::template AnyDense<Derived> &a, Args... args)
    {
        a.setRandom(args...);
    }

    template <typename Derived, typename... Args>
    static void setZero(typename T::template AnyDense<Derived> &a, Args... args)
    {
        a.setZero(args...);
    }
};

using Backend = L_Eigen;

#endif
} // namespace la
using Factory = la::Traits<la::Backend, double>;
using _ = la::Functions<la::Backend, double>;
} // namespace loki