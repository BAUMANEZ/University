//
//  main.cpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#include "../tests/helmholtz/test_helmholtz.hpp"

int main(int argc, const char * argv[]) {
    test_helmholtz::run_omp(false);
    test_helmholtz::run_omp(true);
    return 0;
}
