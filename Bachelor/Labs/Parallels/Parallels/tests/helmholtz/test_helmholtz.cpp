//
//  helmholtz.cpp
//  Parallels
//
//  Created by Арсений Токарев on 14.10.2022.
//

#include "test_helmholtz.hpp"

void test_helmholtz::run(bool paralleled) {
    const size_t t = 50;
    const size_t n =  t*10;
    const double k = double(t)*2.0;
    const double h = 1/double(n - 1);
    const func2d values = [k](double x, double y) -> double {
        return 2*sin(M_PI*y) + pow(k, 2)*(1-x)*x*sin(M_PI*y) + pow(M_PI, 2)*(1-x)*x*sin(M_PI*y);
    };
    const func2d u_00 = [](double x, double y) -> double { return 0.; };
    const func2d u_01 = [](double x, double y) -> double { return 0.; };
    const func2d u_10 = [](double x, double y) -> double { return 0.; };
    const func2d u_11 = [](double x, double y) -> double { return 0.; };

    double start = omp_get_wtime();
    paralleled ? algorithm::omp_helmholtz_red_black(k, h, {0, 1}, {0, 1}, values, {u_00, u_01, u_10, u_11}) : algorithm::helmholtz_red_black(k, h, {0, 1}, {0, 1}, values, {u_00, u_01, u_10, u_11});
    std::cout << "Time=" << omp_get_wtime()-start << "\n\n";

    start = omp_get_wtime();
    paralleled ? algorithm::omp_helmholtz_jacobi(k, h, {0, 1}, {0, 1}, values, {u_00, u_01, u_10, u_11}) : algorithm::helmholtz_jacobi(k, h, {0, 1}, {0, 1}, values, {u_00, u_01, u_10, u_11});
    std::cout << "Time=" << omp_get_wtime()-start << "\n\n";
}
