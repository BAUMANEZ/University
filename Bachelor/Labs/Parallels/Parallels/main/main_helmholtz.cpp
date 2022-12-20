//
//  main.cpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#include "../tests/helmholtz/test_helmholtz.hpp"

int main(int argc, const char * argv[]) {
    test_helmholtz::run(false);
    test_helmholtz::run(true);
    return 0;
}
