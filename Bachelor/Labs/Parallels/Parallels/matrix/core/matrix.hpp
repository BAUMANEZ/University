//
//  matrix.h
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#ifndef matrix_hpp
#define matrix_hpp

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <math.h>
#include "../../helpers/helpers.hpp"

struct matrix final {
private:
    double* values;
    size_t n;
    size_t m;
    
//MARK: - Getters
public:
    inline size_t rows() const { return n; }
    inline size_t cols() const { return m; }

//MARK: - Init
public:
    matrix(size_t n, double initial);
    matrix(size_t n, double l, double r);
    matrix(size_t n, size_t m, double initial);
    matrix(size_t n, size_t m, double l, double r);
    matrix(std::vector<std::vector<double>> nested);
    matrix(const matrix& copy);
    matrix(matrix&& rvalue) noexcept;
    ~matrix();
    
//MARK: - Contents
public:
    double elem(size_t i, size_t j) const;
    double& mutable_elem(size_t i, size_t j);
    
//MARK: - Operators
    matrix operator+(const matrix& rhs);
    matrix operator-(const matrix& rhs);
    matrix operator*(const matrix& rhs);
    void   operator=(const matrix& rhs);
    void   operator=(matrix&& rhs);
    
//MARK: - Linear Algebra Modifications
    void invert();
    matrix inverted() const;
    
//MARK: - Service
public:
    friend std::ostream& operator<<(std::ostream& out, const matrix& ref);
    static void swap(matrix& lhs, matrix& rhs);
    static matrix make_identity(size_t n);
    static matrix make_identity(size_t n, size_t m);
};

#endif /* matrix_h */
