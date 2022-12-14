//
//  norms.cpp
//  Parallels
//
//  Created by Арсений Токарев on 12.10.2022.
//

#include "omp.h"
#include "norms.hpp"

double norms::frobenius(const matrix& lhs, const matrix& rhs) {
    // assert_message(lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols(), "Dimensions are not equal");
    double sum = 0.;
#pragma omp parallel for default(none) shared(lhs, rhs) reduction(+:sum)
    for (size_t i = 0; i < lhs.rows(); ++i)
        for (size_t j = 0; j < lhs.cols(); ++j) {
            const double diff = lhs.elem(i, j)-rhs.elem(i, j);
            sum += diff*diff;
        }
    
    return sqrt(sum);
}
double norms::max_error(const matrix& lhs, const matrix& rhs) {
    // assert_message(lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols(), "Dimensions are not equal");
    double error = 0.;
    for (size_t i = 0; i < lhs.rows(); ++i)
        for (size_t j = 0; j < lhs.cols(); ++j) {
            error = std::max(error, fabs(lhs.elem(i, j)-rhs.elem(i, j)));
        }
    return error;
}
