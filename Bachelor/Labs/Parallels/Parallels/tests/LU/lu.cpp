//
//  lu.cpp
//  Parallels
//
//  Created by Арсений Токарев on 12.10.2022.
//

#include "lu.hpp"

void lu_test::run(method method) {
    matrix test(512, 512, -30., 30.);
    matrix lu_test(test);
    std::cout << "Randomly generated matrix\n";
    std::cout << test;
    clock_t decomposition_time;
    clock_t start = clock();
    switch (method) {
        case non_block:
            std::cout << "LU Decomposition\n";
            algorithm::lu_decomposition(lu_test);
            decomposition_time = clock()-start;
            break;
        case block:
            std::cout << "BLU Decomposition\n";
            algorithm::blu_decomposition(lu_test, 32);
            decomposition_time = clock()-start;
            break;
    }
    std::cout << lu_test;
    std::cout << "Time elapsed on decomposition: " << decomposition_time/double(CLOCKS_PER_SEC) << "\n";
    std::cout << "LU Multiplication\n";
    matrix lu_multiplication(algorithm::lu_multiplication(lu_test));
    std::cout << lu_multiplication;
    std::cout << "Frobenius norm (A, A_{L}*A_{U})=" << norms::frobenius(test, lu_multiplication) << "\n";
    std::cout << "Max error norm (A, A_{L}*A_{U})=" << norms::max_error(test, lu_multiplication) << "\n";
}
