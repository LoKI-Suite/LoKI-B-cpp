//
// Created by daan on 15-5-19.
//

#define EIGEN

#ifndef LOKI_CPP_LINEARALGEBRA_H
#define LOKI_CPP_LINEARALGEBRA_H

#include "Eigen/Dense"

namespace loki {

#ifdef EIGEN

    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::VectorXd Vector;

#endif

}

#endif //LOKI_CPP_LINEARALGEBRA_H
