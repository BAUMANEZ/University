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
            const std::vector<double> i_col = A.col_at(i);
            const std::vector<double> i_row = A.row_at(i);
            for (size_t j = i+1; j < A.cols(); ++j) {
                A.mutable_elem(k, j) -= i_col[k]*i_row[j];
            }
        }
    }
}
void algorithm::blu_decomposition(matrix& A) {
    
}
matrix algorithm::lu_multiplication(const matrix& A) {
    matrix L = matrix::make_identity(A.rows(), A.cols());
    matrix U = matrix::make_identity(A.rows(), A.cols());
    for (size_t i = 0; i < A.rows(); ++i) {
        U.mutable_elem(i, i) = A.elem(i, i);
        for (size_t j = i+1; j < A.cols(); ++j)
            U.mutable_elem(i, j) = A.elem(i, j);
        for (size_t j = 0; j < i; ++j)
            L.mutable_elem(i, j) = A.elem(i, j);
    }
    return L*U;
}
