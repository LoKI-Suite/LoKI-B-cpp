//
// Created by daan on 19-06-2019.
//

#include <cmath>
#include <stdint.h>
#include <cstring>
#include <iostream>

void givens(double a, double b, double &c, double &s) {
    const double root = std::sqrt(a * a + b * b);

    c = a / root;
    s = -b / root;
}

// Matrix A is column major

double *hessenberg(double *A, double *b, uint32_t n) {
    auto *c = new double[n];
    auto *s = new double[n];
    auto *v = new double[n];
    double *x = b;

    std::memcpy(v, A + (n - 1) * n, n * sizeof(double));

    for (uint32_t k = n - 1; k > 0; --k) {
        givens(v[k], A[(k - 1) * n + k], c[k], s[k]);

        x[k] /= c[k] * v[k] - s[k] * A[(k - 1) * n + k];

        for (int64_t i = 0; i < k - 1; ++i) {
            x[i] += x[k] * s[k] * A[(k - 1) * n + i] - x[k] * c[k] * v[i];
            v[i] = c[k] * A[(k - 1) * n + i] + s[k] * v[i];
        }

        x[k - 1] += x[k] * s[k] * A[(k - 1) * n + k - 1] - x[k] * c[k] * v[k - 1];
        v[k - 1] = c[k] * A[(k - 1) * n + k - 1] + s[k] * v[k - 1];
    }


    double t1 = x[0] / v[0];

    for (uint32_t k = 1; k < n; ++k) {
        const double t2 = x[k];
        x[k - 1] = c[k] * t1 - s[k] * t2;
        t1 = c[k] * t2 + s[k] * t1;
    }

    x[n - 1] = t1;

    return x;
}

int main(int argc, char **argv) {
    uint32_t size = 3;

    auto *A = new double[size * size];
    auto *b = new double[size];

    A[0] = 1.3608276348795434e-01;
    A[1] = 1.7341118050852893e+02;
    A[2] = 0.;
    A[3] = 2.3570226039551584e-01;
    A[4] = -4.9394538222198395e+02;
    A[5] = 3.2053408487237232e+02;
    A[6] = 3.0429030972509225e-01;
    A[7] = 3.2053459156859032e+02;
    A[8] = -3.2053459156859032e+02;

    b[0] = 1.;
    b[1] = 0.;
    b[2] = 0.;

    double *x = hessenberg(A, b, size);

    for (uint32_t i = 0; i < size; ++i) {
        std::cout << x[i] << std::endl;
    }

    delete[] A;
    delete[] b;

    return 0;
}