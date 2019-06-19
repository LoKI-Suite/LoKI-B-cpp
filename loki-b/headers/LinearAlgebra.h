//
// Created by daan on 15-5-19.
//

#define EIGEN

#ifndef LOKI_CPP_LINEARALGEBRA_H
#define LOKI_CPP_LINEARALGEBRA_H

#define EIGEN_USE_MKL_ALL

#include "Eigen/Dense"

namespace loki {

#ifdef EIGEN

    typedef Eigen::MatrixXd Matrix;
    // typedef Eigen::internal::BandMatrix<double, Eigen::Dynamic, 1, 1> TridiagonalMatrix;
    typedef Eigen::VectorXd Vector;

    class LinAlg {
    public:
        static void givens(double a, double b, double &c, double &s) {
            const double root = std::sqrt(a * a + b * b);

            c = a / root;
            s = -b / root;
        }

        static double *hessenberg(const double *A, double *b, uint32_t n) {
            auto *c = new double[n];
            auto *s = new double[n];
            auto *v = new double[n];
            double *x = b;

            // Matrix A is column major
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
    };

#endif

}

#endif //LOKI_CPP_LINEARALGEBRA_H
