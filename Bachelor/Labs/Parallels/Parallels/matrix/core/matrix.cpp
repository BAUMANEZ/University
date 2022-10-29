//
//  matrix.cpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#include "matrix.hpp"

//MARK: - Init
matrix::matrix(size_t n, double initial): matrix(n, n, initial) {}
matrix::matrix(size_t n, double l, double r): matrix(n, n, l, r) {}
matrix::matrix(size_t n, size_t m, double initial) {
    this->n = n;
    this->m = m;
    this->values = new double[n*m];
    for (size_t i = 0; i < n*m; ++i) values[i] = initial;
}
matrix::matrix(size_t n, size_t m, double l, double r) {
    this->n = n;
    this->m = m;
    this->values = new double[n*m];
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j) {
            double k = double(rand())/RAND_MAX;
            mutable_elem(i, j) = l+k*(r-l);
        }
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
    this->values = new double[n*m];
    for (size_t i = 0; i < n; ++i) {
        assert_message(nested[i].size() == m, "Number of cols is not equal in each row");
        for (size_t j = 0; j < m; ++j) mutable_elem(i, j) = nested[i][j];
    }
}
matrix::matrix(const matrix& copy) {
    this->n = copy.n;
    this->m = copy.m;
    this->values = new double[n*m];
    std::memcpy(values, copy.values, n*m);
}
matrix::matrix(matrix&& rvalue) noexcept {
    this->n = rvalue.n;
    this->m = rvalue.m;
    delete[] values;
    this->values = rvalue.values;
    rvalue.values = nullptr;
}
matrix::~matrix() {
    delete[] values;
}

//MARK: - Contents
double matrix::elem(size_t i, size_t j) const {
    // assert_message(i<n, "Accessing out of bounds row of matrix. That's illegal");
    // assert_message(j<m, "Accessing out of bounds col of matrix. That's illegal");
    return values[i*m+j];
}
double& matrix::mutable_elem(size_t i, size_t j) {
    // assert_message(i<n, "Accessing out of bounds row of matrix. That's illegal");
    // assert_message(j<m, "Accessing out of bounds col of matrix. That's illegal");
    return values[i*m+j];
}

//MARK: - Operators
matrix matrix::operator+(const matrix& rhs) {
    // assert_message(n == rhs.n && m == rhs.m, "Dimensions are not equal");
    matrix result(*this);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            result.mutable_elem(i, j) += rhs.elem(i, j);
    return result;
}
matrix matrix::operator-(const matrix& rhs) {
    // assert_message(n == rhs.n && m == rhs.m, "Dimensions are not equal");
    matrix result(*this);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            result.mutable_elem(i, j) -= rhs.elem(i, j);
    return result;
}
matrix matrix::operator*(const matrix& rhs) {
    // assert_message(m == rhs.n, "Cols(I) is not equal to Rows(II)");
    matrix result(n, rhs.m, 0.);
    for (size_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < m; ++k) {
            for (size_t j = 0; j < rhs.m; ++j)
                result.mutable_elem(i, j) += elem(i, k)*rhs.elem(k, j);
        }
    }
    return result;
}
void matrix::operator=(const matrix& rhs) {
    n = rhs.n;
    m = rhs.m;
    delete[] values;
    values = new double[n*m];
    std::memcpy(values, rhs.values, n*m);
}
void matrix::operator=(matrix&& rhs) {
    n = rhs.n;
    m = rhs.m;
    delete[] values;
    values = rhs.values;
    rhs.values = nullptr;
}

//MARK: - Linear Algebra Modifications

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
        for (size_t i = 0; i<ref.n; ++i) {
            const double elem = ref.elem(i, j);
            max_j = std::max(max_j, number_of_digits(elem, precision)+((elem < 0. ? 1 : 0)));
        }
        widths[j] = max_j;
    }
    int widths_sum = 0;
    for (auto it = widths.begin(); it < widths.end(); ++it)
        widths_sum += (*it+4);
    
    std::stringstream dashes;
    for (size_t k = 0; k <= widths_sum; ++k) {
        dashes << std::setw(1) << "┅";
    }
    out << dashes.str();
    out << "\n";
    for (size_t i = 0; i < ref.n; ++i) {
        out << "┋ ";
        for (size_t j = 0; j < ref.m; ++j)
            out << std::right << std::setw(widths[j]+1) << cut(ref.elem(i, j), precision) << std::setw(1) << " ┋ ";
        out << "\n";
        out << dashes.str();
        out << "\n";
    }
    out << "\n";
    return out;
}
void matrix::swap(matrix& lhs, matrix& rhs) {
    std::swap(lhs.n, rhs.n);
    std::swap(lhs.m, rhs.m);
    std::swap(lhs.values, rhs.values);
}
matrix matrix::make_identity(size_t n) {
    return matrix::make_identity(n, n);
}
matrix matrix::make_identity(size_t n, size_t m) {
    matrix result(n, m, 0);
    for (size_t i = 0; i < n; ++i)
        result.mutable_elem(i, i) = 1.;
    return result;
}


