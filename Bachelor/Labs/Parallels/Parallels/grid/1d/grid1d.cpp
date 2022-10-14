//
//  grid1d.cpp
//  Parallels
//
//  Created by Арсений Токарев on 14.10.2022.
//

#include "grid1d.hpp"

//MARK: - Init
grid1d::grid1d(double start, double end, double step) {
    this->a = start;
    this->b = end;
    this->h = step;
    this->n = std::max(0, int((end-start)/step) - 1);
}

//MARK: - Operators
double grid1d::operator[](size_t i) const {
    return a+double(i)*h;
}
