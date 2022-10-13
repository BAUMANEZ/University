//
//  algorithm.hpp
//  Parallels
//
//  Created by Арсений Токарев on 08.10.2022.
//

#ifndef algorithm_hpp
#define algorithm_hpp

#include <iostream>
#include "../core/matrix.hpp"

 struct algorithm final {
     static void lu_decomposition(matrix& A);
     static void blu_decomposition(matrix& A, const size_t block);
     static matrix lu_multiplication(const matrix& A);
};

#endif /* algorithm_hpp */
