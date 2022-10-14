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
    #include "./tests/LU/lu.hpp"
//MARK: - when run from terminal/sublime/vscode
#else
    #include "./matrix/core/matrix.cpp"
    #include "./matrix/algorithm/algorithm_lu.cpp"
    #include "./matrix/algorithm/algorithm_helmholtz.cpp"
    #include "./matrix/norms/norms.cpp"
    #include "./tests/LU/lu.cpp"
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
    lu_test::run(lu_test::block, false);
    return 0;
}