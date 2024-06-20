#include "ideas/HessenbergExtra.h"
#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/Log.h"

namespace loki {
namespace LinAlg {

inline bool isUpperHessenberg(const double *A, const uint32_t n)
{
    for (uint32_t k = 0; k < n - 2; ++k)
    {
        for (uint32_t l = k + 2; l < n; ++l)
        {
            if (A[k * n + l] != 0.)
            {
                //std::cerr << l << ", " << k << std::endl;
                return false;
            }
        }
    }

    return true;
}

inline void reduce(double *A, uint32_t first, const uint32_t second, const uint32_t start, const uint32_t n)
{

    if (A[start * n + first] == 0)
        Log<Message>::Error("Could not solve using Hessenberg, please use the LU decomposition.");

    /// \todo Should this result in an exception? Now this is ignored.
    if (second == first)
        std::cerr << "cannot perform reduction on the same row" << std::endl;

    const double ratio = A[start * n + second] / A[start * n + first];

    A[start * n + second] = 0.;
#ifdef LOKIB_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint32_t j = start + 1; j < n; ++j)
    {
        A[j * n + second] -= ratio * A[j * n + first];
    }
}

inline void reduceRowToHess(double *A, const uint32_t a, const uint32_t b, const uint32_t start, const uint32_t n)
{
    for (uint32_t i = a; i < b; ++i)
    {
        reduce(A, i, b, start + i - a, n);
    }
}

/* TODO: [Optional] Update Hessenberg reduction such that it also alters the b vector.
 *
 * Idea: Use the 'bandgap' between superelastic diagonals to optimize the algorithm.
 *
 * 1. c[i] - c[i-1] determines the maximum number of subtractions to be saved by subtracting
 * any lower row from c[i].
 *
 * Now let k > i
 *
 * 2. c[k] - c[k-1] determines the maximum number of subtractions that can be saved by
 * subtracting with c[k].
 * 3. c[k+1] - c[k] determines the number of rows j, for c[i]<=j<c[i+1], following c[i]
 * that subtraction with c[k]+j can be applied to.
 *
 * Thus for each c[i] we would like to find the most optimal row to subtract. Thus we
 * want to choose a block c[k] -- c[k+1] to be subtracted that skips the maximum number
 * of subtractions possible. This value can be calculated by multiplying the number of
 * subtractions that we skip per subtraction with the number of subtractions that can
 * be applied, hence:
 *
 * The optimality factor 'of' can be calculated through:
 *
 * (min(c[i] - c[i-1], c[k] - c[k-1]) - 1) * min(c[i+1] - c[i], c[k+1] - c[k])
 *
 * or
 *
 * (min(g[i-1], g[k-1]) - 1) * min(g[i], g[k])
 */

inline void hessenbergReductionOptimal(double *A, const uint32_t *c, uint32_t n, uint32_t cn)
{
    // Array storing the gaps between two subsequent bands in band array (c),
    // such that g[i] stores c[i+1] - c[i].
    std::vector<uint32_t> gv(cn-1);
    auto *g = gv.data();

    for (uint32_t k = 0; k < cn - 1; ++k)
    {
        g[k] = c[k + 1] - c[k];
    }

    // The element c[0] should be set to 1 coinciding with the first subdiagonal of the matrix.
    // Therefore the loop should start at i = 1.

    for (uint32_t i = 1; i < cn - 1; ++i)
    {
        // g[i-1] is at least 1, in which case no substractions can be saved.
        if (g[i - 1] == 1)
        {
            // Row reduce the block c[i] through c[i+1].
            for (uint32_t j = 0; j < g[i]; ++j)
            {
                //                        std::cerr << "first hess reduce" << std::endl;
                reduceRowToHess(A, j + 1, j + c[i], j, n);
            }
        }
        else
        {
            uint32_t ofMax = 0, of = 0, kOpt = 0;
            uint32_t optRows = 0, optSavesPerRow = 0;

            for (uint32_t k = i + 1; k < cn; ++k)
            {
                const uint32_t numRows = std::min(g[i], g[k]);
                const uint32_t numSavesPerRow = std::min(g[i - 1], g[k - 1]);

                of = numRows * numSavesPerRow;

                if (of > ofMax)
                {
                    ofMax = of;
                    kOpt = k;
                    optRows = numRows;
                    optSavesPerRow = numSavesPerRow;
                }
            }

            for (uint32_t j = 0; j < optRows; ++j)
            {
                // Smart reduce
                //                        std::cerr << "smart reduce" << std::endl;
                reduce(A, c[kOpt] + j, c[i] + j, j, n);

                // If the row is not yet Hessenberg, then reduce.
                // In other words, if j + optSavesPerRow (the position of the element that we
                // have eliminated + the number of substractions that we can skip) is smaller
                // than the index of the element on the first subdiagonal (c[i]-1) then the
                // row is not yet upper Hessenberg and we need to reduce it further (in a
                // top-down approach).
                if (optSavesPerRow < c[i] - 1)
                {
                    //                            std::cerr << "hess reduce after smart reduce" << std::endl;
                    reduceRowToHess(A, j + optSavesPerRow + 1, j + c[i], j + optSavesPerRow, n);
                }
            }

            // TODO: CHECK BOUNDS
            for (uint32_t j = optRows; j < g[i]; ++j)
            {
                //                        std::cerr << "optRows: " << optRows << std::endl;
                reduceRowToHess(A, j + 1, j + c[i], j, n);
            }
        }
    }

    // This part is correct, since we cannot reduce the rows below the last
    // band in a bottom up fashion.
    const uint32_t lastRow = c[cn - 1];

    //            std::cerr << "last loop" << std::endl;
    for (uint32_t i = lastRow; i < n; ++i)
    {
        reduceRowToHess(A, i - lastRow + 1, i, i - lastRow, n);
    }
}

inline void hessenbergReduction(double *A, const uint32_t *c, uint32_t n, uint32_t cn)
{
    for (uint32_t ci = 0; ci < cn - 1; ++ci)
    {
        for (uint32_t row = c[ci]; row < c[ci + 1]; ++row)
        {
            reduceRowToHess(A, row - c[ci] + 1, row, row - c[ci], n);
        }
    }

    const uint32_t lastRow = c[cn - 1];

    for (uint32_t i = lastRow; i < n; ++i)
    {
        reduceRowToHess(A, i - lastRow + 1, i, i - lastRow, n);
    }
}

inline void swapRows(uint32_t *p, const uint32_t one, const uint32_t two)
{
    uint32_t temp = p[one];

    p[one] = p[two];
    p[two] = temp;
}

inline void reduceRowToHess(double *A, uint32_t *p, const uint32_t a, const uint32_t b, const uint32_t start,
                            const uint32_t n)
{
    // a is the row at with the first element at the subdiagonal
    // b is the row that we want to reduce
    // thus row a is used to reduce row b

    for (uint32_t i = a; i < b; ++i)
    {
        // Before we reduce a row, we find the row starting that has the highest value starting element
        // (thus the row needs to have its first element at position 'col', underneath the sub diagonal)
        // and place this row such that its first element is on the sub diagonal
        const uint32_t col = start + i - a;

        double pivMax = 0, piv = 0;
        uint32_t iMax = 0;

        // start searching at the subdiagonal and continue until the row that we want to reduce

        for (uint32_t j = col + 1; j <= b; ++j)
        {
            piv = std::abs(A[col * n + p[j]]);

            if (piv > pivMax)
            {
                pivMax = piv;
                iMax = j;
            }
        }

        if (iMax != col + 1)
        {
            swapRows(p, iMax, col + 1);
            std::cerr << "swapping\n" << std::endl;
        }

        reduce(A, p[i], p[b], col, n);
    }
}

void hessenbergReductionPartialPiv(double *A, const uint32_t *c, uint32_t *p, uint32_t n, uint32_t cn)
{
    // p is the permutation vector (it keeps track of the current positions of the rows. Thus p[i] is the current
    // position of the row that was originally at position i.
    for (uint32_t ci = 0; ci < cn - 1; ++ci)
    {
        for (uint32_t row = c[ci]; row < c[ci + 1]; ++row)
        {

            reduceRowToHess(A, p, row - c[ci] + 1, row, row - c[ci], n);
        }
    }

    const uint32_t lastRow = c[cn - 1];

    for (uint32_t i = lastRow; i < n; ++i)
    {
        reduceRowToHess(A, p, i - lastRow + 1, i, i - lastRow, n);
    }
}

double *hessenberg(double *A, double *b, const uint32_t *p, uint32_t n, HessenbergWorkspace& hws)
{
    if (hws.size()<n)
    {
        Log<Message>::Error("hessenberg: workspace is too small.");
    }
    auto *c = hws.c.data();
    auto *s = hws.s.data();
    auto *v = hws.v.data();
    double *x = b;

    std::memcpy(v, A + (n - 1) * n, n * sizeof(double));

    for (uint32_t i = n - 1; i > 0; --i)
    {
        const uint32_t k = p[i];

        impl::givens(v[k], A[(k - 1) * n + k], c[k], s[k]);

        x[k] /= c[k] * v[k] - s[k] * A[(k - 1) * n + k];

        for (int64_t j = 0; j < k - 1; ++j)
        {
            x[j] += x[k] * s[k] * A[(k - 1) * n + j] - x[k] * c[k] * v[j];
            v[j] = c[k] * A[(k - 1) * n + j] + s[k] * v[j];
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

double *hessenberg(double *A, double *b, const uint32_t *p, uint32_t n)
{
    HessenbergWorkspace hws(n);
    return hessenberg(A,b,p,n,hws);
}

} // namespace LinAlg
} // namespace loki
