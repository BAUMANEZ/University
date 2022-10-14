//
//  algorithm.hpp
//  Parallels
//
//  Created by Арсений Токарев on 08.10.2022.
//

#ifndef algorithm_hpp
#define algorithm_hpp

#include <iostream>
#include <array>
#include <functional>
//#include "omp.h"
#include "../core/matrix.hpp"
#include "../norms/norms.hpp"
#include "../../grid/1d/grid1d.hpp"
#include "../../helpers/helpers.hpp"

 struct algorithm final {
//MARK: - LU
     static void lu_decomposition(matrix& A);
     static void omp_lu_decomposition(matrix& A);
     static void blu_decomposition(matrix& A, const size_t block);
     static void omp_blu_decomposition(matrix& A, const size_t block);
     static matrix lu_multiplication(const matrix& A);
     
//MARK: - Helmholtz
 public:
     static void helmholtz(double k, std::array<double, 2> abscissa, std::array<double, 2> ordinate, func2d values, std::array<func2d, 4> boundary);
};

#endif /* algorithm_hpp */
