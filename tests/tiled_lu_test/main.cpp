#include <iostream>

#include <chrono>

#include "TiledMatrix.h"
#include "LinearAlgebra.h"

using namespace la;

typedef Eigen::VectorXd Vector;

int main(int argc, char **argv)
{
    uint32_t block_size, dim;

    if (argc != 3)
    {
        // std::cerr << "./tiled-lu <block-size> <problem-dimension>" << std::endl;
        // return -1;
        block_size = 128;
        dim = 1024;
    }
    else
    {
        block_size = std::stoi(std::string(argv[1]));
        dim = std::stoi(std::string(argv[2]));
    }

    Eigen::MatrixXd matrix(dim, dim), c_matrix;

    Vector b = Vector::Zero(dim);
    b[0] = 1;

    matrix.setRandom();

    c_matrix = matrix;

    TiledMatrix<uint32_t, double> t_matrix(block_size);

    auto begin = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();

    for (uint32_t i = 0; i < 10; ++i)
    {
        t_matrix.fromDense(matrix);

        begin = std::chrono::high_resolution_clock::now();

        tiled_lu(t_matrix);

        end = std::chrono::high_resolution_clock::now();
        std::cerr << "Tiled:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;

        begin = std::chrono::high_resolution_clock::now();

        matrix.partialPivLu();

        end = std::chrono::high_resolution_clock::now();
        std::cerr << "Eigen:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;
    }
}
