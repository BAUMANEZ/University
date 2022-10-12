//
//  lu.cpp
//  Parallels
//
//  Created by Арсений Токарев on 12.10.2022.
//

#include "lu.hpp"

void lu_test::run(method method, bool parallel) {
    matrix test(4, 4, -5., 5);
    matrix lu_test(test);
    std::cout << "Randomly generated matrix\n";
    std::cout << test;
    switch (method) {
        case non_block:
            std::cout << "LU Decomposition\n";
            algorithm::lu_decomposition(lu_test);
            break;
        case block:
            std::cout << "BLU Decomposition\n";
            algorithm::blu_decomposition(lu_test);
            break;
    }
    std::cout << lu_test;
    matrix lu_multiplication(algorithm::lu_multiplication(lu_test));
    std::cout << "LU Multiplication\n";
    std::cout << lu_multiplication;
    std::cout << "Frobenius norm (A, A_{L}*A_{U})=" << norms::frobenius(test, lu_multiplication) << "\n";
    std::cout << "Max error norm (A, A_{L}*A_{U})=" << norms::max_error(test, lu_multiplication) << "\n";
}
