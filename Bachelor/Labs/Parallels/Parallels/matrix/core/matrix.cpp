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
    (this->values).resize(n*m, initial);
}
matrix::matrix(size_t n, size_t m, double l, double r) {
    this->n = n;
    this->m = m;
    (this->values).resize(n*m);
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
    (this->values).reserve(n*m);
    for (auto it = nested.begin(); it<nested.end(); ++it) {
        assert_message((*it).size()==m, "Number of cols is not equal in each row");
        (this->values).insert((this->values.end()), (*it).begin(), ((*it).end()));
    }
}
matrix::matrix(const matrix& copy) {
    this->n = copy.n;
    this->m = copy.m;
    this->values = copy.values;
}
matrix::matrix(matrix&& rvalue) noexcept {
    this->n = rvalue.n;
    this->m = rvalue.m;
    this->values = std::move(rvalue.values);
}

//MARK: - Setters
void matrix::resize(size_t n, double fill) {
    resize(n, n, fill);
}
void matrix::resize(size_t n, size_t m, double fill) {
    const matrix copy(*this);
    this->n = n;
    this->m = m;
    values = std::vector<double>(n*m, fill);
    const size_t min_n = std::min(n, copy.n);
    const size_t min_m = std::min(m, copy.m);
    for (size_t i = 0; i < min_n; ++i) {
        const size_t copy_row = i*copy.m;
        const size_t new_row = i*m;
        std::copy(
            copy.values.begin()+copy_row,
            copy.values.begin()+copy_row+min_m,
            values.begin()+new_row
        );
    }
}

//MARK: - Contents
double matrix::elem(size_t i, size_t j) const {
    assert_message(i<n, "Accessing out of bounds row of matrix. That's illegal");
    assert_message(j<m, "Accessing out of bounds col of matrix. That's illegal");
    return values[i*m+j];
}
double& matrix::mutable_elem(size_t i, size_t j) {
    assert_message(i<n, "Accessing out of bounds row of matrix. That's illegal");
    assert_message(j<m, "Accessing out of bounds col of matrix. That's illegal");
    return values[i*m+j];
}
std::vector<double> matrix::row_at(size_t i) const {
    std::vector<double> result;
    std::copy(values.begin()+i*m, values.begin()+(i+1)*m, std::back_inserter(result));
    return result;
}
std::vector<double> matrix::col_at(size_t j) const {
    std::vector<double> result;
    result.reserve(m);
    for (auto it = values.begin(); it < values.end(); it += m)
        result.push_back(*(it+j));
    return result;
}
matrix matrix::submatrix(size_t n, size_t m, size_t i, size_t j) const {
    assert_message(i+n<=this->n, "Submatrix has more rows than the original");
    assert_message(j+m<=this->m, "Submatrix has more cols than the original");
    matrix result(n, m, 0.);
    for (size_t k = i; k < i+n; ++k)
        for (size_t l = j; l < j+m; ++l)
            result.mutable_elem(k-i, l-j) = elem(k, l);
    return result;
}

//MARK: Operators
matrix matrix::operator+(const matrix& rhs) {
    assert_message(n == rhs.n && m == rhs.m, "Dimensions are not equal");
    matrix result(*this);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            result.mutable_elem(i, j) += rhs.elem(i, j);
    return result;
}
matrix matrix::operator-(const matrix& rhs) {
    assert_message(n == rhs.n && m == rhs.m, "Dimensions are not equal");
    matrix result(*this);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
            result.mutable_elem(i, j) -= rhs.elem(i, j);
    return result;
}
matrix matrix::operator*(const matrix& rhs) {
    assert_message(m == rhs.n, "Cols(I) is not equal to Rows(II)");
    matrix result(n, rhs.m, 0.);
    for (size_t i = 0; i < n; ++i) {
//        const std::vector<double> i_row = row_at(i);
        for (size_t k = 0; k < m; ++k) {
//            const std::vector<double> i_col = rhs.row_at(k);
            for (size_t j = 0; j < rhs.m; ++j)
                result.mutable_elem(i, j) += elem(i, k)*rhs.elem(k, j);//i_row[k]*i_col[j];
        }
    }
    return result;
}
void matrix::operator=(const matrix& rhs) {
    n = rhs.n;
    m = rhs.m;
    values = rhs.values;
}
void matrix::operator=(const matrix&& rhs) {
    n = rhs.n;
    m = rhs.m;
    values = std::move(rhs.values);
}

//MARK: - Linear Algebra Modifications
void invert() {
    
}
//matrix inverted();


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
        for (size_t j = 0; j < ref.m; ++j) {
            std::stringstream stream;
            stream << ref.elem(i, j);
            const std::string string = stream.str();
            const size_t endIndex = string.size() <= precision ? string.size() : precision;
            out << std::right << std::setw(widths[j]+1) << string.substr(0, endIndex) << std::setw(1) << " ┋ ";
        }
        out << "\n";
        out << dashes.str();
        out << "\n";
    }
    out << "\n";
    return out;
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


