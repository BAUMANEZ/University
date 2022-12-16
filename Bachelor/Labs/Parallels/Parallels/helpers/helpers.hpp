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

std::string cut(double number, size_t n);

bool is_odd(size_t n);

int number_of_digits(double num, int precision = 16);

#endif /* helpers_hpp */
