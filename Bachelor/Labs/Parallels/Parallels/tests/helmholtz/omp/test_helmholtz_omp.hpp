//
//  helmholtz.hpp
//  Parallels
//
//  Created by Арсений Токарев on 14.10.2022.
//

#ifndef helmholtz_omp_hpp
#define helmholtz_omp_hpp

#include <iostream>
#include <math.h>
#include "../../../helpers/helpers.hpp"
#include "../../../matrix/algorithm/algorithm.hpp"

struct test_helmholtz_omp final {
    static void run(bool paralleled);
};

#endif /* helmholtz_hpp */
