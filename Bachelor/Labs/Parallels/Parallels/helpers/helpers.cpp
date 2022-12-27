//
//  helpers.cpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#include "helpers.hpp"

std::string cut(double number, size_t n) {
    std::stringstream stream;
    stream << number;
    const std::string string = stream.str();
    const size_t endIndex = string.size() <= n ? string.size() : n;
    return string.substr(0, endIndex);
}

bool is_odd(size_t n) {
   return (n & 1);
}

int number_of_digits(double num, int precision) {
    int digits = 0;
    double original = num;

    //before point
    long num2=num;
    while(num2>0)
    {
        digits++;
        num2=num2/10;
    }
    if(original==0)
        digits=1;
    num = original;
    double no_float;
    no_float = original*(pow(10, (16-digits)));

    //after point
    long long int total=(long long int)no_float;
    int no_of_digits, extrazeroes=0;
    for(int i=0; i<16; i++)
    {
        int dig;
        dig=total%10;
        total=total/10;
        if(dig!=0)
            break;
        else
            extrazeroes++;
    }
    no_of_digits=16-extrazeroes;
    return std::min(no_of_digits, precision);
}

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

bool is_first_or_last(int index, int size) {
    return (index == 0) || (index == size - 1);
}
