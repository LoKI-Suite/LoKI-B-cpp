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

#endif

}

#endif //LOKI_CPP_LINEARALGEBRA_H
