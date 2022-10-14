//
//  algorithm.hpp
//  Parallels
//
//  Created by Арсений Токарев on 08.10.2022.
//

#ifndef algorithm_hpp
#define algorithm_hpp

#include <iostream>
#include "omp.h"
#include "../core/matrix.hpp"

 struct algorithm final {
//MARK: - LU
     static void lu_decomposition(matrix& A);
     static void omp_lu_decomposition(matrix& A);
     static void blu_decomposition(matrix& A, const size_t block);
     static void omp_blu_decomposition(matrix& A, const size_t block);
     static matrix lu_multiplication(const matrix& A);
     
//MARK: - Helmholtz
 public:
     static void helmholtz();
};

#endif /* algorithm_hpp */
