//
// Created by daan on 30-4-19.
//

#include <iostream>
#include <vector>

#define EIGEN

class Matrix;

#ifdef EIGEN
class EigenMatrix {

};

typedef std::vector<double> Vector;

class Matrix : public EigenMatrix {
public:
    Vector solve(const Vector &b) {
        return {};
    }

    void shout() {
        std::cout << "EigenMatrix" << std::endl;
    }
};
#endif

#ifdef ARMADILLO
class ArmadilloMatrix {

};

typedef std::vector<double> Vector;

class Matrix : public ArmadilloMatrix {
public:
    Vector solve(const Vector &b) {
        return {};
    }

    void shout() {
        std::cout << "ArmadilloMatrix" << std::endl;
    }
};
#endif

int main (int argc, char ** argv)
{
    Vector test_vector;

    Matrix matrix;

    matrix.shout();

    return 0;
}
