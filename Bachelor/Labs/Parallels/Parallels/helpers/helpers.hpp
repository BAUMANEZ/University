//
//  helpers.hpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#ifndef helpers_hpp
#define helpers_hpp

#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <functional>
#include "math.h"

#define assert_message(exp, msg) assert(((void)msg, exp))

typedef std::function<double(double, double)> func2d;

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

inline int number_of_digits(double num, int precision = 16) {
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

#endif /* helpers_hpp */
