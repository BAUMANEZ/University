//
//  main.cpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#include <iostream>

//MARK: - when run from xcode
#if DEBUG || RELEASE
    #include "./matrix/core/matrix.hpp"
    #include "./matrix/algorithm/algorithm.hpp"
    #include "./matrix/norms/norms.hpp"
    #include "./tests/lu/test_lu.hpp"
    #include "./tests/helmholtz/test_helmholtz.hpp"
//MARK: - when run from terminal/sublime/vscode
#else
    #include "omp.h"
    #include "./matrix/core/matrix.cpp"
    #include "./grid/1d/grid1d.cpp"
    #include "./matrix/algorithm/algorithm_lu.cpp"
    #include "./matrix/algorithm/algorithm_helmholtz.cpp"
    #include "./matrix/norms/norms.cpp"
    #include "./tests/lu/test_lu.cpp"
    #include "./tests/helmholtz/test_helmholtz.cpp"
#endif

std::ostream& operator<<(std::ostream& out, const std::vector<double>& ref) {
    const int precision = 8;
    out << "Vector, len=(" << ref.size() << "), adress=" << &ref << "\n";
    if (!(ref.size() < 10)) {
        return out;
    }
    int width = 1;
    for (const auto& x: ref)
        width = std::max(width, number_of_digits(x, precision));
    for (const auto& x: ref)
        out << std::left << std::setw(2) << "|" << std::right << std::setw(width) << x << std::setw(3) <<"|\n";
    out << "\n";
    return out;
}

int main(int argc, const char * argv[]) {
    test_lu::run(test_lu::non_block, false);
    test_lu::run(test_lu::non_block, true);
    test_lu::run(test_lu::block, false);
    test_lu::run(test_lu::block, true);

    // test_helmholtz::run(false);
    // test_helmholtz::run(true);
    return 0;
}
