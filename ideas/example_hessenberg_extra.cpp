#include "ideas/HessenbergExtra.h"
#include "LoKI-B/LinearAlgebra.h"
#include <iostream>
#include <stdexcept>

using loki::Matrix;
using loki::Vector;

void testHessenberg(const Matrix& Aorig)
{
    if (Aorig.rows()!=Aorig.cols())
    {
        throw std::runtime_error("matrix is not square.");
    }
    Matrix A(Aorig);
    A.row(0).setZero();
    A(0,0)=A(1,1);
    const Matrix::Index dim = A.rows();
    Vector b(dim); 
    b.setZero();
    b[0] = 1.;

    std::cout << "A:\n" << A << std::endl;
    std::cout << "b:\n" << b.transpose() << std::endl;
    std::cout << "Ainv:\n" << A.inverse() << std::endl;
    std::cout << "Ainv*b:\n" << (A.inverse()*b).transpose() << std::endl;
    // 'standard' Hessenberg.

    Vector x1(b);
    loki::LinAlg::hessenberg(A.data(), x1.data(), dim);
    std::cout << "x1:    " << x1.transpose() << std::endl;

#if 0
    /** \todo Below are some snippets that were previously in ElectronKinetics.cpp,
     *  updated so they compile cleanly. If we want to keep the permuting hessenberg
     *  implementation, we should complete and enable this code.
     */

    /** \todo Is it correct that p is not used after the call to
     *  hessenbergReductionPartialPiv?
     */
    std::vector<uint32_t> p(dim);
    for (Matrix::Index i = 0; i != dim; ++i)
    {
        p[i] = i;
    }

    /// \todo What should this be? The diagonals containing non-zero entries?
    std::vector<uint32_t> superElasticThresholds;
    // note: keeping this empty makes the code crash
    superElasticThresholds.push_back(0);
    
    loki::LinAlg::hessenbergReductionPartialPiv(A.data(), superElasticThresholds.data(), p.data(),
      dim, superElasticThresholds.size());
   
    Vector x2(b);
    loki::LinAlg::hessenberg(A.data(), x2.data(), p.data(), dim);

    std::cout << "x2:    " << x2.transpose() << std::endl;
    std::cout << "x1-x2: " << (x1-x2).transpose() << std::endl;
#endif
}

int main()
{
   loki::Matrix A(10,10);
   A.diagonal().fill(1.0);
   A.diagonal(-1).fill(-1.0);
   testHessenberg(A);
}
