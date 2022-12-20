//
//  main_helmoltz_mpi.cpp
//  Parallels
//
//  Created by Арсений Токарев on 20.12.2022.
//

#include "../tests/helmholtz/test_helmholtz.hpp"

int main(int argc, const char * argv[]) {
    test_helmholtz::run_mpi(false);
    test_helmholtz::run_mpi(true);
    return 0;
}
