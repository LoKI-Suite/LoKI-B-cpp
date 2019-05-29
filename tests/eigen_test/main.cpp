//
// Created by daan on 8-5-19.
//

#include <iostream>
#include <Eigen/Dense>

#define EIGEN

namespace loki {
    class eigen;

    template<typename Library>
    class Traits;

    template<>
    class Traits<eigen> {
    public:
        typedef Eigen::MatrixXd matrix_type;
        typedef Eigen::VectorXd c_vec_type;

        // This expression might be wrong. Also, this approach might not work at all
        // (due to the nature of the expression templates).
        typedef const Eigen::Solve<Eigen::LLT<matrix_type>, c_vec_type> solve_type;
    };

    template<typename Library>
    class T_Matrix {
    protected:
        typename Traits<Library>::matrix_type *matrix{};

    public:
        // pure virtual solve function with solve type available in trait class
        virtual typename Traits<Library>::solve_type solve(const typename Traits<Library>::c_vec_type &b) = 0;

        virtual void fill() = 0;

        T_Matrix() = default;

        virtual ~T_Matrix() { delete matrix; }
    };

    /*
     * There are now two options, either we specialize the template class for each library. Or
     * we keep T_Matrix as an abstract class and for each library (in an #ifdef), we declare
     * a Matrix class that inherits T_Matrix<libraryname> and implement the functions there.
     */

#ifdef EIGEN

    typedef Traits<eigen> LinAlg;

    class Matrix : public T_Matrix<eigen> {
    public:
        explicit Matrix(uint32_t dim) {
            this->matrix = new Traits<eigen>::matrix_type(dim, dim);
        }

        // This function does not work because of expression templates.
        // I tried to fix it by defining a 'solve_type' in the trait class that should hold
        // the resulting expression template, but it does not function correctly.
        // One could make it work by passing the output vector by reference as an argument.
        // However I am curious how to solve this problem in a more elegant way.
        Traits<eigen>::solve_type solve(const Traits<eigen>::c_vec_type &b) override {

            return this->matrix->llt().solve(b);
        }

        void fill() override {
            this->matrix->setRandom();
        }

        ~Matrix() override = default;
    };

#endif

    /*
     * We can also immediately inherit from our Trait class matrix type. However,
     * the user now still has access to the library specific functions which is
     * probably not desirable.
     */

    template<typename Library>
    class B_Matrix : public Traits<Library>::matrix_type {
    public:
        template<typename... Args>
        explicit B_Matrix(Args... args) : Traits<Library>::matrix_type(args...) {}

        virtual ~B_Matrix() = default;

        virtual typename Traits<Library>::solve_type solve(const typename Traits<Library>::c_vec_type &b) = 0;

        virtual void fill() = 0;
    };

    class Alt_Matrix : public B_Matrix<eigen> {
    public:
        explicit Alt_Matrix(uint32_t dim) : B_Matrix(dim, dim) {}

        ~Alt_Matrix() override = default;

        Traits<eigen>::solve_type solve(const Traits<eigen>::c_vec_type &b) override {
            return this->llt().solve(b);
        }

        void fill() override {
            this->setRandom();
        }
    };
}

/*
 * Another method to show that it also does not work when using the auto keyword.
 */
auto testFunc(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) -> decltype(A.llt().solve(b)) {
    return A.llt().solve(b);
}

auto addFunc(const Eigen::VectorXd &a, const Eigen::VectorXd &b) {
    return a + b;
}

template <typename derivedA, typename derivedB>
const Eigen::Solve<Eigen::LLT<derivedA, 1>, derivedB>
solve(const Eigen::MatrixBase<derivedA> &A, const Eigen::MatrixBase<derivedB> &b) {
    return A.llt().solve(b);
}

template <typename derivedA, typename derivedB>
void solve(const Eigen::MatrixBase<derivedA> &A, Eigen::MatrixBase<derivedB> &x, const Eigen::MatrixBase<derivedB> &b) {
    x = A.llt().solve(b);
}

int main(int argc, char **argv) {
    Eigen::MatrixXd A(10, 10);
    A.fill(1);
    Eigen::VectorXd b(10);
    b.fill(1);
    Eigen::VectorXd x(10);

    solve(A, x, b);

    // All three implementations will cause std::bad_alloc().
//    x = testFunc(A, b);

//    loki::Matrix matrix(10);
//    matrix.fill();
//    x = matrix.solve(b);

//    loki::Alt_Matrix altMatrix(10);
//    altMatrix.fill();
//    x = altMatrix.solve(b);

    // However this works:
//    x = addFunc(b, b);


    std::cout << x << std::endl;

    return 0;
}