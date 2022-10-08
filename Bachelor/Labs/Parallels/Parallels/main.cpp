//
//  main.cpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#include <iostream>
#include "matrix/matrix.hpp"

int main(int argc, const char * argv[]) {
    matrix a({
        {1.333, 2.934782, 3.333333},
        {4.235235, 5.3214, 6.323232},
        {7.124312, 8.3242, 9.3328389}
    });
    
    std::cout << a;
    return 0;
}
