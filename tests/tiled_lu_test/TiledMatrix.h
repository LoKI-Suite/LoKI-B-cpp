#ifndef TILED_MATRIX_H
#define TILED_MATRIX_H

#define EIGEN_USE_LAPACKE
#define EIGEN_USE_BLAS
// #define EIGEN_USE_MKL_ALL

#include <Eigen/Dense>
#include <vector>
#include <cassert>

namespace la
{
template <typename IndexType, typename ElementType>
class TiledMatrix
{
    typedef Eigen::Matrix<ElementType, Eigen::Dynamic, Eigen::Dynamic> Dense;

public:
    TiledMatrix(IndexType block_size) : block_size(block_size)
    {
    }

    Dense &get_tile(const IndexType i, const IndexType j)
    {
        assert(i < num_tile_col && j < num_tile_row);

        return tiles[j * num_tile_col + i];
    }

    void fromDense(const Dense &matrix)
    {
        num_row = matrix.rows();
        num_col = matrix.cols();

        num_tile_row = (IndexType)std::ceil(matrix.rows() / (float)block_size);
        num_tile_col = (IndexType)std::ceil(matrix.cols() / (float)block_size);

        const IndexType num_elem_col_full_blocks = (num_tile_col - 1) * block_size;
        const IndexType leftover_col_dim = matrix.cols() - num_elem_col_full_blocks;

        const IndexType num_elem_row_full_blocks = (num_tile_row - 1) * block_size;
        const IndexType leftover_row_dim = matrix.rows() - num_elem_row_full_blocks;

        tiles.resize(num_tile_row * num_tile_col);

#pragma omp parallel for
        for (uint32_t tile_col = 0; tile_col < num_tile_col - 1; ++tile_col)
        {
            for (uint32_t tile_row = 0; tile_row < num_tile_row - 1; ++tile_row)
            {
                get_tile(tile_row, tile_col) = matrix.block(tile_row * block_size, tile_col * block_size, block_size, block_size);
            }

            get_tile(num_tile_row - 1, tile_col) = matrix.block(num_elem_row_full_blocks, tile_col * block_size, leftover_row_dim, block_size);
        }

#pragma omp parallel for
        for (uint32_t tile_row = 0; tile_row < num_tile_row - 1; ++tile_row)
        {
            get_tile(tile_row, num_tile_col - 1) = matrix.block(tile_row * block_size, num_elem_col_full_blocks, block_size, leftover_col_dim);
        }

        get_tile(num_tile_row - 1, num_tile_col - 1) = matrix.block(num_elem_row_full_blocks, num_elem_col_full_blocks, leftover_row_dim, leftover_col_dim);
    }

    IndexType tile_rows()
    {
        return num_tile_row;
    }

    IndexType tile_cols()
    {
        return num_tile_col;
    }

    IndexType rows()
    {
        return num_row;
    }

    IndexType cols()
    {
        return num_col;
    }

    IndexType get_block_size()
    {
        return block_size;
    }

private:
    IndexType block_size;
    IndexType num_tile_col, num_tile_row, num_col, num_row;

    std::vector<Dense, Eigen::aligned_allocator<Dense>> tiles;
};
} // namespace la

#endif