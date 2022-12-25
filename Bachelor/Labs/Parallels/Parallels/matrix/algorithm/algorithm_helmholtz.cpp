//
//  algorithm_helmholtz.cpp
//  Parallels
//
//  Created by Арсений Токарев on 13.10.2022.
//

#include "omp.h"
#include "algorithm.hpp"

#define eps 1e-7

size_t algorithm::helmholtz_red_black(
    matrix& result,
    const double k, 
    const double h, 
    const grid1d& x,
    const grid1d& y,
    const func2d& values
) {
    const double div = 1./(4.+pow(k,2)*pow(h, 2));
    const double gamma = pow(h, 2)*div;
    size_t iterations = 0;
    matrix previous(result);
    do {
        ++iterations;
        matrix::swap(result, previous);
#pragma omp parallel for default(none) shared(result, previous, gamma, div, x, y, values)
        for (size_t j = 1; j < y.steps()-1; ++j)
            for (size_t i = is_odd(j) ? 1 : 2; i < x.steps()-1; i += 2)
                 result.mutable_elem(j, i) = gamma*values(x[i], y[j])+div*(previous.elem(j+1, i)+previous.elem(j-1, i)+previous.elem(j, i+1)+previous.elem(j, i-1));
#pragma omp parallel for default(none) shared(result, gamma, div, x, y, values)
        for (size_t j = 1; j < y.steps()-1; ++j)
            for (size_t i = is_odd(j) ? 2 : 1; i < x.steps()-1; i += 2)
                 result.mutable_elem(j, i) = gamma*values(x[i], y[j])+div*(result.elem(j+1, i)+result.elem(j-1, i)+result.elem(j, i+1)+result.elem(j, i-1));
    } while (norms::frobenius(result, previous) > eps);

    return iterations;
}

size_t algorithm::helmholtz_jacobi(
    matrix& result,
    const double k, 
    const double h, 
    const grid1d& x,
    const grid1d& y,
    const func2d& values
) {
    const double div = 1./(4.+pow(k,2)*pow(h, 2));
    const double gamma = pow(h, 2)*div; 
    size_t iterations = 0;
    matrix previous(result);
    do {
        ++iterations;
        matrix::swap(result, previous);
#pragma omp parallel for default(none) shared(result, previous, gamma, div, x, y, values)
        for (size_t j = 1; j < y.steps()-1; ++j)
            for (size_t i = 1; i < x.steps()-1; ++i)
                result.mutable_elem(j, i) = gamma*values(x[i], y[j])+div*(previous.elem(j+1, i)+previous.elem(j-1, i)+previous.elem(j, i+1)+previous.elem(j, i-1));
    } while (norms::frobenius(result, previous) > eps);

    return iterations;
}
