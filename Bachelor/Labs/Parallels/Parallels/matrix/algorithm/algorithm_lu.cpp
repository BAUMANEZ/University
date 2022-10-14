//
//  algorithm.cpp
//  Parallels
//
//  Created by Арсений Токарев on 08.10.2022.
//

#include "algorithm.hpp"

void algorithm::lu_decomposition(matrix& A) {
    size_t min_dim = std::min(A.rows(), A.cols());
    for (size_t i = 0; i < min_dim; ++i) {
        //MARK: (1)
        const double div = 1./A.elem(i, i);
        for (size_t k = i+1; k < A.rows(); ++k)
            A.mutable_elem(k, i) *= div;
        
        //MARK: (2)
        for (size_t k = i+1; k < A.rows(); ++k) {
//            const std::vector<double> i_col = A.col_at(i);
//            const std::vector<double> i_row = A.row_at(i);
            for (size_t j = i+1; j < A.cols(); ++j) {
                A.mutable_elem(k, j) -= A.elem(k, i)*A.elem(i, j);//i_col[k]*i_row[j];
            }
        }
    }
}
void algorithm::omp_lu_decomposition(matrix& A) {
    size_t min_dim = std::min(A.rows(), A.cols());
    for (size_t i = 0; i < min_dim; ++i) {
        const double div = 1./A.elem(i, i);
#pragma omp parallel for default(none) shared(A, i, div)
            for (size_t k = i+1; k < A.rows(); ++k)
                A.mutable_elem(k, i) *= div;
#pragma omp parallel for default(none) shared(A, i)
        for (size_t k = i+1; k < A.rows(); ++k) {
            for (size_t j = i+1; j < A.cols(); ++j) {
                A.mutable_elem(k, j) -= A.elem(k, i)*A.elem(i, j);
            }
        }
    }
}
void algorithm::blu_decomposition(matrix& A, const size_t block) {
    assert_message(A.rows() == A.cols(), "Matrix must be squared");
    matrix submatrix(A.rows(), block, 0.);
    for (size_t i = 0; i < A.rows()-1; i += block) {
        //MARK: (1)
        /// - Result is an LU block of size [(i, i+n), (i, i+block)]
        /// - L22 part is [(i, i+block), (i, i+block)]
        /// - L32 part is [(i+block, n), (i, i+block)]
        for (size_t k = i; k < A.rows(); ++k)
            for (size_t l = i; l < i+block; ++l)
                submatrix.mutable_elem(k-i, l-i) = A.elem(k, l);
        algorithm::lu_decomposition(submatrix);
        for (size_t k = i; k < A.rows(); ++k)
            for (size_t l = i; l < i+block; ++l)
                A.mutable_elem(k, l) = submatrix.elem(k-i, l-i);
        
        //MARK: (2)
        /// - U23 part is [(i, i+block), (i+block, n)]
        /// - U23 is modified by multiplying Inversed(L22) by original U23 =>
        /// - In fact L22 is a lower triangular matrix with 1's on the diaginal =>
        /// - Therefore k goes from i+1 (i_th row is equal to original U23)
        for (size_t k = i+1; k < i+block; ++k)
            for (size_t l = i+block; l < A.rows(); ++l)
                for (size_t m = i; m < k; ++m)
                    A.mutable_elem(k, l) -= A.elem(k, m)*A.elem(m, l); //U23[k, l] -= L22[k, m]*U23[m, l]
        
        //MARK: (3)
        /// - A33 part is [(i+block, n), (i+block, n)]
        for (size_t k = i+block; k < A.rows(); ++k)
            for (size_t m = i; m < i+block; ++m)
                for (size_t l = i+block; l < A.rows(); ++l)
                    A.mutable_elem(k, l) -= A.elem(k, m)*A.elem(m, l); //A33[k, l] -= L32[k, m]*U23[m, l]
    }
}
void algorithm::omp_blu_decomposition(matrix& A, const size_t block) {
    assert_message(A.rows() == A.cols(), "Matrix must be squared");
    matrix submatrix(A.rows(), block, 0.);
    for (size_t i = 0; i < A.rows()-1; i += block) {
        for (size_t k = i; k < A.rows(); ++k)
            for (size_t l = i; l < i+block; ++l)
                submatrix.mutable_elem(k-i, l-i) = A.elem(k, l);
        algorithm::lu_decomposition(submatrix);
        for (size_t k = i; k < A.rows(); ++k)
            for (size_t l = i; l < i+block; ++l)
                A.mutable_elem(k, l) = submatrix.elem(k-i, l-i);
        
        for (size_t k = i+1; k < i+block; ++k)
            for (size_t l = i+block; l < A.rows(); ++l)
                for (size_t m = i; m < k; ++m)
                    A.mutable_elem(k, l) -= A.elem(k, m)*A.elem(m, l);
        
#pragma omp parallel for default(none) shared(A, i, block)
        for (size_t k = i+block; k < A.rows(); ++k)
            for (size_t m = i; m < i+block; ++m)
                for (size_t l = i+block; l < A.rows(); ++l)
                    A.mutable_elem(k, l) -= A.elem(k, m)*A.elem(m, l);
    }
}
matrix algorithm::lu_multiplication(const matrix& A) {
    matrix L = matrix::make_identity(A.rows(), A.cols());
    matrix U = matrix::make_identity(A.rows(), A.cols());
    for (size_t i = 0; i < A.rows(); ++i) {
        U.mutable_elem(i, i) = A.elem(i, i);
        for (size_t j = 0; j < i; ++j)
            L.mutable_elem(i, j) = A.elem(i, j);
        for (size_t j = i+1; j < A.cols(); ++j)
            U.mutable_elem(i, j) = A.elem(i, j);
    }
    return L*U;
}
