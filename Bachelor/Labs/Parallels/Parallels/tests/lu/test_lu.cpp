//
//  lu.cpp
//  Parallels
//
//  Created by Арсений Токарев on 12.10.2022.
//

#include "omp.h"
#include "test_lu.hpp"

void test_lu::run(method method, bool paralleled) {
    omp_set_num_threads(paralleled ? 18 : 1);
    std::cout << (paralleled ? "PARALLELED" : "SEQUENTIAL") << "\n\n";
    for (size_t n = 512; n <= 8192; n *= 2) {
        matrix lu_test(n, n, -1000., 1000.);
        // matrix lu_test(test);
        std::cout << "Random ";
        std::cout << lu_test;
        switch (method) {
            case non_block:
                double start = omp_get_wtime();
                algorithm::lu_decomposition(lu_test);
                const double end = omp_get_wtime()-start;
                std::cout << "LU Decomposition: time=" << end << "\n";
                break;
            case block:
                std::cout << "BLU Decomposition: ";
                for (size_t block = 32; block <= 32; block += 32) {
                    double start = omp_get_wtime();
                    algorithm::blu_decomposition(lu_test, block);
                    const double end = omp_get_wtime()-start;
                    std::cout << "block=" << block << ", time=" << end << ";  ";
                }
                std::cout << "\n";
                break;
        }
        // if (n == 512) {
        //     std::cout << "LU Multiplication: ";
        //     matrix lu_multiplication(algorithm::lu_multiplication(lu_test));
        //     std::cout << lu_multiplication;
        //     std::cout << "Frobenius norm (A, A_{L}*A_{U})=" << norms::frobenius(test, lu_multiplication) << "\n";
        //     std::cout << "Max error norm (A, A_{L}*A_{U})=" << norms::max_error(test, lu_multiplication) << "\n";
        // }
        std::cout << "\n\n";
    }
    std::cout << "END OF " << (paralleled ? "PARALLELED" : "SEQUENTIAL") << "\n\n";
}
