//
//  algorithm_helmholtz.cpp
//  Parallels
//
//  Created by Арсений Токарев on 13.10.2022.
//

#include "algorithm.hpp"

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
        for (size_t i = 1; i < y.steps()-1; ++i)
            for (size_t j = is_odd(i) ? 1 : 2; j < x.steps()-1; j += 2)
                 result.mutable_elem(i, j) = gamma*values(x[i], y[j])+div*(previous.elem(i+1, j)+previous.elem(i-1, j)+previous.elem(i, j+1)+previous.elem(i, j-1));
#pragma omp parallel for default(none) shared(result, gamma, div, x, y, values)
        for (size_t i = 1; i < y.steps()-1; ++i)
            for (size_t j = is_odd(i) ? 2 : 1; j < x.steps()-1; j += 2)
                 result.mutable_elem(i, j) = gamma*values(x[i], y[j])+div*(result.elem(i+1, j)+result.elem(i-1, j)+result.elem(i, j+1)+result.elem(i, j-1));
    } while (norms::frobenius(result, previous) > __DBL_EPSILON__);

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
        for (size_t i = 1; i < x.steps()-1; ++i)
            for (size_t j = 1; j < y.steps()-1; ++j)
                result.mutable_elem(i, j) = gamma*values(x[i], y[j])+div*(previous.elem(i+1, j)+previous.elem(i-1, j)+previous.elem(i, j+1)+previous.elem(i, j-1));
    } while (norms::frobenius(result, previous) > __DBL_EPSILON__);

    return iterations;
}
