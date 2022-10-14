//
//  algorithm_helmholtz.cpp
//  Parallels
//
//  Created by Арсений Токарев on 13.10.2022.
//

#include "algorithm.hpp"

void algorithm::helmholtz(
    double k,
    std::array<double, 2> abscissa,
    std::array<double, 2> ordinate,
    func2d values,
    std::array<func2d, 4> boundary
) {
    const double h = 0.1;
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
    
    //MARK: Iteration Method
    const double div = 1./(4.+pow(k,2)*pow(h, 2));
    const double gamma = pow(h, 2)/div;
    matrix previous(result);
    size_t iterations = 0;
    do {
        ++iterations;
        previous = result;
        for (size_t i = 1; i < x.steps()-1; ++i) {
            for (size_t j = 1; j < y.steps()-1; ++j) {
                result.mutable_elem(i, j) = gamma*values(x[i], y[j]) + div*(previous.elem(i+1, j)+previous.elem(i-1, j)+previous.elem(i, j+1)+previous.elem(i, j-1));
            }
        }
    } while (norms::frobenius(previous, result) > __DBL_EPSILON__);
    std::cout << "Number of iterations=" << std::max(0, (int)iterations-1) << "\n";
}
