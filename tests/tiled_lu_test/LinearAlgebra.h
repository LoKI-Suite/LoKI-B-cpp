#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include "TiledMatrix.h"

namespace la
{
// Assumes a square matrix A of size n*n
template <typename IndexType, typename ElementType>
static void lu_decomposition_parallel(ElementType *A, const IndexType n)
{
    for (IndexType index = 0; index < n - 1; ++index)
    {
#pragma omp parallel for
        for (IndexType col = index + 1; col < n; ++col)
        {
            A[col * n + index] /= A[index * n + index];

            for (IndexType row = index + 1; row < n; ++row)
                A[col * n + row] -= A[index * n + row] * A[col * n + index];
        }
    }
}

// Assumes a square tiled matrix
template <typename IndexType, typename ElementType>
void tiled_lu(TiledMatrix<IndexType, ElementType> &matrix)
{
    typedef Eigen::Matrix<ElementType, Eigen::Dynamic, Eigen::Dynamic> Dense;

    const IndexType n = matrix.tile_rows();
    const IndexType bs = matrix.get_block_size();

    for (IndexType index = 0; index < n - 1; ++index)
    {
        lu_decomposition_parallel(matrix.get_tile(index, index).data(), bs);

#pragma omp parallel for
        for (IndexType col = index + 1; col < n; ++col)
        {
            matrix.get_tile(index, col) = matrix.get_tile(index, index).template triangularView<Eigen::Lower>().solve(matrix.get_tile(index, col));
        }

#pragma omp parallel for
        for (IndexType row = index + 1; row < n; ++row)
        {
            matrix.get_tile(row, index) = matrix.get_tile(index, index).template triangularView<Eigen::UnitUpper>().template solve<Eigen::OnTheRight>(matrix.get_tile(row, index));
        }

#pragma omp parallel for collapse(2)
        for (uint32_t col = index + 1; col < n; ++col)
        {
            for (uint32_t row = index + 1; row < n; ++row)
                matrix.get_tile(row, col) -= matrix.get_tile(row, index) * matrix.get_tile(index, col);
        }
    }

    lu_decomposition_parallel(matrix.get_tile(n - 1, n - 1).data(), matrix.rows() - (n - 1) * bs);
}
} // namespace la

#endif