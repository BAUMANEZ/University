//
//  helmholtz.hpp
//  Parallels
//
//  Created by Арсений Токарев on 14.10.2022.
//

#ifndef helmholtz_hpp
#define helmholtz_hpp

#include <iostream>
#include <math.h>
#include "../../helpers/helpers.hpp"
#include "../../matrix/algorithm/algorithm.hpp"

struct test_helmholtz final {
    static void run_omp(bool paralleled);
    static void run_mpi(bool paralleled);
};

#endif /* helmholtz_hpp */
