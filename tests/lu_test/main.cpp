//
// Created by daan on 6-6-19.
//

#include <iostream>

void LUDecomposition(const double *A, double *L, double *U, const uint32_t dim);

int main (int argc, char ** argv)
{
    uint32_t dim = 2;

    auto * A = new double[dim*dim]{3, 1, 4, 2};
    auto * U = new double[dim*dim];
    auto * L = new double[dim*dim];

    LUDecomposition(A, L, U, dim);

    for (uint32_t i = 0; i < dim; ++i) {
        for (uint32_t j = 0; j < dim; ++j) {
            std::cout << U[i * dim + j] << '\t';
        }

        std::cout << '\n';
    }
}

void LUDecomposition(const double *A, double *L, double *U, const uint32_t dim) {
    for (uint32_t k = 0; k < dim; ++k) {
        for (uint32_t m = 0; m < dim; ++m) {

            uint32_t index = k * dim + m;

            U[index] = A[index];

            for (uint32_t j = 0; j < k; ++j) {
                U[index] -= L[k * dim + j]*U[j * dim + m];
            }
        }

        const double Ukk = U[k * dim + k];

        for (uint32_t i = 0; i < dim; ++i) {
            if (i == k) {
                L[i * dim + i] = 1;
                continue;
            }

            const uint32_t index = i * dim + k;

            L[index] = A[index];

            for (uint32_t j = 0; j < k; ++j) {
                L[index] -= L[i * dim + j] * U[j * dim + k];
            }

            L[index] /= Ukk;
        }
    }
}