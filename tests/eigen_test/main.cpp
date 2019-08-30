//
// Created by daan on 8-5-19.
//

#include <iostream>
#include "backend_abstraction.h"

#define EIGEN

using namespace loki;

int main(int argc, char **argv)
{
    Factory::Matrix A;
    _::setRandom(A, 10, 10);

    Factory::Vector b;
    _::setRandom(b, 10);

    Factory::Vector x(10);

    x = _::pLUSolve(A, b);

    std::cout << x << std::endl;

    return 0;
}