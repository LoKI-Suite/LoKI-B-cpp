//
// Created by daan on 24-06-2019.
//

#include <cstdint>
#include <iostream>

void print(const double *A, const uint32_t *p, uint32_t n);

void reduce(double *A, uint32_t first, const uint32_t second, const uint32_t start,
            const uint32_t n) {

    if (A[start * n + first] == 0)
        std::cerr << "Could not solve using Hessenberg, please use the LU decomposition.\n";

    if (second == first)
        std::cerr << "cannot perform reduction on the same row" << std::endl;

    const double ratio = A[start * n + second] / A[start * n + first];

    A[start * n + second] = 0.;

    for (uint32_t j = start + 1; j < n; ++j) {
        A[j * n + second] -= ratio * A[j * n + first];
    }
}

void reduceRowToHess(double *A, const uint32_t a, const uint32_t b, const uint32_t start,
                     const uint32_t n) {
    for (uint32_t i = a; i < b; ++i) {
        reduce(A, i, b, start + i - a, n);
    }
}

void swapRows(uint32_t *p, const uint32_t one, const uint32_t two) {
    uint32_t temp = p[one];

    p[one] = p[two];
    p[two] = temp;
}

void reduceRowToHess(double *A, uint32_t *p, const uint32_t a, const uint32_t b, const uint32_t start,
                     const uint32_t n) {
    // a is the row at with the first element at the subdiagonal
    // b is the row that we want to reduce\
    // thus row a is used to reduce row b

    for (uint32_t i = a; i < b; ++i) {
        // Before we reduce a row, we have to find the row with the same starting element (underneath the sub diagonal)
        // and place this row such that its first element is on the sub diagonal
        const uint32_t col = start + i - a;


        double pivMax = 0, piv = 0;
        uint32_t iMax = 0;

        // start searching at the subdiagonal and continue until the row that we want to reduce

        for (uint32_t j = col + 1; j <= b; ++j) {
            piv = std::abs(A[col * n + p[j]]);

            if (piv > pivMax) {
                pivMax = piv;
                iMax = j;
            }
        }

        if (iMax != col + 1)
            swapRows(p, iMax, col + 1);

        print(A, p, n);

        reduce(A, p[i], p[b], col, n);

        std::cout << '\n';
        print(A, p, n);
        std::cout << '\n';
    }
}

void hessenbergReduction(double *A, const uint32_t *c, uint32_t n, uint32_t cn) {
    for (uint32_t ci = 0; ci < cn - 1; ++ci) {
        for (uint32_t row = c[ci]; row < c[ci + 1]; ++row) {
            reduceRowToHess(A, row - c[ci] + 1, row, row - c[ci], n);
        }
    }

    const uint32_t lastRow = c[cn - 1];

    for (uint32_t i = lastRow; i < n; ++i) {
        reduceRowToHess(A, i - lastRow + 1, i, i - lastRow, n);
    }
}

void hessenbergReductionPartialPiv(double *A, const uint32_t *c, uint32_t *p, uint32_t n, uint32_t cn) {
    // p is the permutation vector (it keeps track of the current positions of the rows. Thus p[i] is the current
    // position of the row that was originally at position i.
    for (uint32_t ci = 0; ci < cn - 1; ++ci) {
        for (uint32_t row = c[ci]; row < c[ci + 1]; ++row) {

            reduceRowToHess(A, p, row - c[ci] + 1, row, row - c[ci], n);
        }
    }

    const uint32_t lastRow = c[cn - 1];

    for (uint32_t i = lastRow; i < n; ++i) {
        reduceRowToHess(A, p, i - lastRow + 1, i, i - lastRow, n);
    }
}

void print(const double *A, const uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j < n; ++j) {
            std::cout << A[j * n + i] << '\t';
        }

        std::cout << std::endl;
    }
}

void print(const double *A, const uint32_t *p, const uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j < n; ++j) {
            std::cout << A[j * n + p[i]] << '\t';
        }

        std::cout << std::endl;
    }
}

int main(int argc, char **argv) {

    const uint32_t n = 6;

    auto *A = new double[n * n]{0.};

    A[0] = 1.;
    A[1] = 1.;
    A[7] = 1.;
    A[8] = 1.;
    A[14] = 1.;
    A[15] = 1.;
    A[21] = 1.;
    A[22] = 1.;
    A[28] = 1.;
    A[29] = 1.;
    A[35] = 1.;

    A[32] = 2.;

    // One band at the third subdiagonal.
    auto *c = new uint32_t[3]{3, 4, 5};

    A[3] = 3.;
    A[10] = 10.;
    A[17] = 1.;

    A[4] = 8.;
    A[11] = 1.;

    A[5] = 1.;

    // vector to keep track of the row-wise pivots
    auto *p = new uint32_t[n];

    for (uint32_t k = 0; k < n; ++k)
        p[k] = k;

    print(A, n);

    hessRedOther(A, c, n, 3);

    std::cout << "\nAfter reduction\n";

    print(A, p, n);

    std::cout << "\nPermutation vector\n";

    for (uint32_t i = 0; i < n; ++i) {
        std::cout << p[i] << std::endl;
    }

    return 0;
}