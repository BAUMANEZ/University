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
        //(1)
        const double div = 1./A.elem(i, i);
        for (size_t k = i+1; k < A.rows(); ++k)
            A.mutable_elem(k, i) *= div;
        
        //(2)
        for (size_t k = i+1; k < A.rows(); ++k) {
//            const std::vector<double> i_col = A.col_at(i);
//            const std::vector<double> i_row = A.row_at(i);
            for (size_t j = i+1; j < A.cols(); ++j) {
                A.mutable_elem(k, j) -= A.elem(k, i)*A.elem(i, j);//i_col[k]*i_row[j];
            }
        }
    }
}
void algorithm::blu_decomposition(matrix& A, const size_t block) {
    assert_message(A.rows() == A.cols(), "Matrix must be squared");
    matrix L22(block, block, 0.);
    matrix L32(A.rows()-block, block, 0.);
    matrix U23(block, A.rows()-block, 0.);
    matrix submatrix(A.rows(), block, 0.);
    for (size_t i = 0; i < A.rows()-1; i+= block) {
        //(1)
        for (size_t k = i; k < A.rows(); ++k)
            for (size_t l = i; l < i+block; ++l)
                submatrix.mutable_elem(k-i, l-i) = A.elem(k, l);
        algorithm::lu_decomposition(submatrix);
        for (size_t k = i; k < A.rows(); ++k)
            for (size_t l = i; l < i+block; ++l)
                A.mutable_elem(k, l) = submatrix.elem(k-i, l-i);
        
        //(2.1)
        for (size_t k = i; k < i+block; ++k) {
            L22.mutable_elem(k-i, k-i) = 1.;
            for (size_t l = i; l < k; ++l)
                L22.mutable_elem(k-i, l-i) = A.elem(k, l);
        }
        for (size_t k = i+block; k < A.rows(); ++k) {
            for (size_t l = i; l < i+block; ++l)
                L32.mutable_elem(k-(i+block), l-i) = A.elem(k, l);
        }
        
        //(2.2)
        for (size_t k = i; k < i+block; ++k) {
            for (size_t l = i+block; l < A.rows(); ++l)
                U23.mutable_elem(k-i, l-(i+block)) = A.elem(k, l);
        }
        for (size_t k = 1; k < block; ++k) {
            for (size_t l = 0; l < A.rows()-(i+block); ++l)
                for (size_t m = 0; m < k; ++m)
                    U23.mutable_elem(k, l) -= L22.elem(k, m)*U23.elem(m, l);
        }
        for (size_t k = i; k < i+block; ++k) {
            for (size_t l = i+block; l < A.rows(); ++l)
                A.mutable_elem(k, l) = U23.elem(k-i, l-(i+block));
        }
        
        //(3)
        for (size_t k = i+block; k < A.rows(); ++k) {
            for (size_t m = 0; m < block; ++m) {
                for (size_t l = i+block; l < A.rows(); ++l)
                    A.mutable_elem(k, l) -= L32.elem(k-(i+block), m)*U23.elem(m, l-(i+block));
            }
        } 
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
