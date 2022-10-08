//
//  matrix.cpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#ifndef matrix_cpp
#define matrix_cpp

#include "matrix.hpp"
#include <math.h>

//MARK: - Getters
size_t matrix::rows() const {
    return this->n;
}
size_t matrix::cols() const {
    return this->m;
}

//MARK: - Init
matrix::matrix(size_t n) {
    matrix(n, n);
}
matrix::matrix(size_t n, double initial) {
    matrix(n, n, initial);
}
matrix::matrix(size_t n, size_t m) {
    this->n = n;
    this->m = m;
    values.reserve(n*m);
}
matrix::matrix(size_t n, size_t m, double initial) {
    this->n = n;
    this->m = m;
    values = std::vector<double>(n*m, initial);
}
matrix::matrix(std::vector<std::vector<double>> nested) {
    if (nested.empty()) {
        this->values = {};
        this->n = 0;
        this->m = 0;
        return;
    }
    this->n = nested.size();
    this->m = nested[0].size();
    (this->values).reserve(n*m);
    for (auto it = nested.begin(); it<nested.end(); ++it) {
        assert_message((*it).size()==m, "Number of cols is not equal in each row");
        (this->values).insert((this->values.end()), (*it).begin(), ((*it).end()));
    }
}

//MARK: - Contents
double matrix::elem(size_t i, size_t j) const {
    assert_message(i<n, "Accessing out of bounds row of matrix. That's illegal");
    assert_message(j<m, "Accessing out of bounds col of matrix. That's illegal");
    return values[i*n+j];
}
double& matrix::mutable_elem(size_t i, size_t j) {
    assert_message(i<n, "Accessing out of bounds row of matrix. That's illegal");
    assert_message(j<m, "Accessing out of bounds col of matrix. That's illegal");
    return values[i*n+j];
}

//MARK: - Service
std::ostream& operator<<(std::ostream& out, const matrix& ref) {
    const int precision = 8;
    out << "MATRIX, dim=(" << ref.n << ", " << ref.m << "), adress=" << &ref << "\n";
    if (!(ref.n < precision*3) || !(ref.m < precision*3)) {
        return out;
    }
    std::vector<int> widths(ref.m, 1);
    for (size_t j = 0; j<ref.m; ++j) {
        int max_j = 1;
        for (size_t i = 0; i<ref.n; ++i)
            max_j = std::max(max_j, number_of_digits(ref.elem(i, j), precision));
        widths[j] = max_j;
    }
    int widths_sum = 3;
    for (auto it = widths.begin(); it < widths.end(); ++it)
        widths_sum += (*it+3);
    
    std::string dashes = " ";
    for (size_t k = 0; k <= widths_sum; ++k) {
        dashes += "┅";
    }
    out << dashes;
    out << "\n";
    for (size_t i = 0; i < ref.n; ++i) {
        out << " ┋ ";
        for (size_t j = 0; j < ref.m; ++j)
            out << std::setprecision(precision) << std::right
                << std::setw(widths[j]+1) << ref.elem(i, j) << " ┋ ";
        out << "\n";
        out << dashes;
        out << "\n";
    }
    out << "\n";
    return out;
}


#endif
