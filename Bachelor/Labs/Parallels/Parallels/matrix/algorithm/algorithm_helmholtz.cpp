//
//  algorithm_helmholtz.cpp
//  Parallels
//
//  Created by Арсений Токарев on 13.10.2022.
//

#include "algorithm.hpp"

void algorithm::helmholtz_red_black(
    const double k,
    const double h,
    std::array<double, 2> abscissa,
    std::array<double, 2> ordinate,
    func2d values,
    std::array<func2d, 4> boundary
) {
    const grid1d x(abscissa[0], abscissa[1], h);
    const grid1d y(ordinate[0], ordinate[1], h);
    matrix result(x.steps(), y.steps(), 0.);
    for (size_t i = 0; i < x.steps(); ++i) {
        result.mutable_elem(i, 0) = boundary[0](x[i], y[0]);
        result.mutable_elem(i, y.steps()-1) = boundary[1](x[i], y[y.steps()-1]);
    }
    for (size_t j = 0; j < y.steps(); ++j) {
        result.mutable_elem(0, j) = boundary[2](x[0], y[j]);
        result.mutable_elem(x.steps()-1, j) = boundary[3](x[x.steps()-1], y[j]);
    }

    matrix previous(result);
    matrix test(x.steps(), y.steps(), 0.);
    for (size_t i = 1; i < x.steps()-1; ++i)
        for (size_t j = 1; j < y.steps()-1; ++j)
            test.mutable_elem(i, j) = (1-x[i])*x[i]*sin(M_PI*y[j]);

    const double div = 1./(4.+pow(k,2)*pow(h, 2));
    const double gamma = pow(h, 2)*div;
    size_t iterations = 0;
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
    
    std::cout << "Frobenius norm(Red-black)=" << norms::frobenius(result, test) << "\n";
    std::cout << "Max error norm(Red-black)=" << norms::max_error(result, test) << "\n";
    std::cout << "Number of iterations(Red-black)=" << std::max(0, (int)iterations-1) << "\n";
}
void algorithm::helmholtz_jacobi(
    const double k,
    const double h,
    std::array<double, 2> abscissa,
    std::array<double, 2> ordinate,
    func2d values,
    std::array<func2d, 4> boundary
) {
    const grid1d x(abscissa[0], abscissa[1], h);
    const grid1d y(ordinate[0], ordinate[1], h);
    matrix result(x.steps(), y.steps(), 0.);
    for (size_t i = 0; i < x.steps(); ++i) {
        result.mutable_elem(i, 0) = boundary[0](x[i], y[0]);
        result.mutable_elem(i, y.steps()-1) = boundary[1](x[i], y[y.steps()-1]);
    }
    for (size_t j = 0; j < y.steps(); ++j) {
        result.mutable_elem(0, j) = boundary[2](x[0], y[j]);
        result.mutable_elem(x.steps()-1, j) = boundary[3](x[x.steps()-1], y[j]);
    }

    matrix previous(result);
    matrix test(x.steps(), y.steps(), 0.);
    for (size_t i = 1; i < x.steps()-1; ++i)
        for (size_t j = 1; j < y.steps()-1; ++j)
            test.mutable_elem(i, j) = (1-x[i])*x[i]*sin(M_PI*y[j]);
    const double div = 1./(4.+pow(k,2)*pow(h, 2));
    const double gamma = pow(h, 2)*div;
    
    size_t iterations = 0;
    do {
        ++iterations;
        matrix::swap(result, previous);
#pragma omp parallel for default(none) shared(result, previous, gamma, div, x, y, values)
        for (size_t i = 1; i < x.steps()-1; ++i)
            for (size_t j = 1; j < y.steps()-1; ++j)
                result.mutable_elem(i, j) = gamma*values(x[i], y[j])+div*(previous.elem(i+1, j)+previous.elem(i-1, j)+previous.elem(i, j+1)+previous.elem(i, j-1));
    } while (norms::frobenius(result, previous) > __DBL_EPSILON__);
    std::cout << "Frobenius norm(Jacobi)=" << norms::frobenius(result, test) << "\n";
    std::cout << "Max error norm(Red-black)=" << norms::max_error(result, test) << "\n";
    std::cout << "Number of iterations(Jacobi)=" << std::max(0, (int)iterations-1) << "\n";
}
