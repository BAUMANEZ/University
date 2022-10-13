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
    std::vector<double> values;
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
    
//MARK: - Setters
    void resize(size_t n, double fill);
    void resize(size_t n, size_t m, double fill);
    
//MARK: - Contents
public:
    double elem(size_t i, size_t j) const;
    double& mutable_elem(size_t i, size_t j);
    std::vector<double> row_at(size_t i) const;
    std::vector<double> col_at(size_t i) const;
    matrix submatrix(size_t n, size_t m, size_t i = 0, size_t j = 0) const;
    
//MARK: Operators
    matrix operator+(const matrix& rhs);
    matrix operator-(const matrix& rhs);
    matrix operator*(const matrix& rhs);
    void   operator=(const matrix& rhs);
    void   operator=(const matrix&& rhs);
    
//MARK: - Linear Algebra Modifications
    void invert();
    matrix inverted() const;
    
//MARK: - Service
public:
    friend std::ostream& operator<<(std::ostream& out, const matrix& ref);
    static matrix make_identity(size_t n);
    static matrix make_identity(size_t n, size_t m);
};

#endif /* matrix_h */
