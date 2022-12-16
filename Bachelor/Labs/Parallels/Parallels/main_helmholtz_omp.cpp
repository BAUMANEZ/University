//
//  main.cpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#include <iostream>
#include "./tests/helmholtz/test_helmholtz.hpp"

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
    test_helmholtz::run(false);
    test_helmholtz::run(true);
    return 0;
}
