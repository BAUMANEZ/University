//
//  matrix.h
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#ifndef matrix_hpp
#define matrix_hpp

#include <stdarg.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "../helpers/helpers.hpp"

struct matrix {
private:
    std::vector<double> values;
    size_t n;
    size_t m;
    
//MARK: - Getters
public:
    size_t rows() const;
    size_t cols() const;
    
//MARK: - Setters
    void set_size(size_t n);
    void set_size(size_t n, size_t m);

//MARK: - Init
public:
    matrix(size_t n);
    matrix(size_t n, double initial);
    matrix(size_t n, size_t m);
    matrix(size_t n, size_t m, double initial);
    matrix(std::vector<std::vector<double>> nested);
    
//MARK: - Contents
public:
    double elem(size_t i, size_t j) const;
    double& mutable_elem(size_t i, size_t j);
    
//MARK: - Service
public:
    friend std::ostream& operator<<(std::ostream& out, const matrix& ref);
    static matrix identity_matrix(size_t n);
    static matrix identity_matrix(size_t n, size_t m);
};

#endif /* matrix_h */
