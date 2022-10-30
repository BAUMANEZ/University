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
#include "omp.h"
#include "../core/matrix.hpp"
#include "../norms/norms.hpp"
#include "../../grid/1d/grid1d.hpp"
#include "../../helpers/helpers.hpp"

 struct algorithm final {
//MARK: - LU
    static void lu_decomposition(matrix& A);
    static void blu_decomposition(matrix& A, const size_t block);
    static matrix lu_multiplication(const matrix& A);
     
//MARK: - Helmholtz
 public:
    static size_t helmholtz_red_black(
        matrix& result,
        const double k, 
        const double h, 
        const grid1d& x,
        const grid1d& y,
        const func2d& values
    );
   static size_t helmholtz_jacobi(
        matrix& result,
        const double k, 
        const double h, 
        const grid1d& x,
        const grid1d& y,
        const func2d& values
    );
};

#endif /* algorithm_hpp */
