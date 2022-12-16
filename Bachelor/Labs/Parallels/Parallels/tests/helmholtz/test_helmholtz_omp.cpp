//
//  helmholtz.cpp
//  Parallels
//
//  Created by Арсений Токарев on 14.10.2022.
//

#include <omp.h>
#include "test_helmholtz.hpp"

void test_helmholtz::run(bool paralleled) {
    std::cout << (paralleled ? "PARALLELED" : "SEQUENTIAL") << "\n\n";
    omp_set_num_threads(paralleled ? 18 : 1);

    for (size_t i = 1; i <= 2; ++i) {
        const size_t t = 320*i;
        const size_t n =  t*10;
        const double h = 1./double(n - 1);
        const double k = double(t)*20.0;
        const func2d values = [k](double x, double y) -> double {
            return 2*sin(M_PI*y) + pow(k, 2)*(1-x)*x*sin(M_PI*y) + pow(M_PI, 2)*(1-x)*x*sin(M_PI*y);
        };
        const func2d u_00 = [](double x, double y) -> double { return 0.; };
        const func2d u_01 = [](double x, double y) -> double { return 0.; };
        const func2d u_10 = [](double x, double y) -> double { return 0.; };
        const func2d u_11 = [](double x, double y) -> double { return 0.; };
        const grid1d x(0., 1., h);
        const grid1d y(0., 1., h);
        const std::array<func2d, 4> boundary = {u_00, u_01, u_10, u_11};

        matrix test(x.steps(), y.steps(), 0.);
        for (size_t i = 0; i < x.steps(); ++i) {
            test.mutable_elem(i, 0) = boundary[0](x[i], y[0]);
            test.mutable_elem(i, y.steps()-1) = boundary[1](x[i], y[y.steps()-1]);
        }
        for (size_t j = 0; j < y.steps(); ++j) {
            test.mutable_elem(0, j) = boundary[2](x[0], y[j]);
            test.mutable_elem(x.steps()-1, j) = boundary[3](x[x.steps()-1], y[j]);
        }
        matrix analytical(x.steps(), y.steps(), 0.);
        for (size_t i = 1; i < x.steps()-1; ++i)
            for (size_t j = 1; j < y.steps()-1; ++j)
                analytical.mutable_elem(i, j) = (1-x[i])*x[i]*sin(M_PI*y[j]);

        matrix copy_test(test);

        std::cout << "===========================";
        std::cout << "\n";

        std::cout << "N=" << n << ", h=" << h << ", k=" << k << "\n\n";

        //MARK: Red-Black
        std::cout << "***RED-BLACK***:\n";
        double start = omp_get_wtime();
        size_t iterations_red_black = algorithm::helmholtz_red_black(test, k, h, x, y, values);
        std::cout << "Time=" << omp_get_wtime()-start << "\n";
        std::cout << "Iterations=" << std::max(0, (int)iterations_red_black-1) << ", ||frobenius||=" << norms::frobenius(test, analytical);


        std::cout << "\n\n";

        //MARK: Jacobi
        std::cout << "***JACOBI***:\n";
        start = omp_get_wtime();
        size_t iterations_jacobi = algorithm::helmholtz_jacobi(copy_test, k, h, x, y, values);
        std::cout << "Time=" << omp_get_wtime()-start << "\n";
        std::cout << "Iterations=" << std::max(0, (int)iterations_jacobi-1) << ", ||frobenius||=" << norms::frobenius(copy_test, analytical);

        std::cout << "\n";
        std::cout << "===========================";
        std::cout << "\n\n";
    }

    std::cout << "END OF " << (paralleled ? "PARALLELED" : "SEQUENTIAL") << "\n";
    std::cout << "\n\n";
}
