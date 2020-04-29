#include <cmath>
#include "Utils.hpp"

std::vector<double> add_vectors(std::vector<double>& x, std::vector<double>& y, double scaling_factor) {
    if (x.size() != y.size()) {
        throw std::length_error("Utils.add_vectors: x and y have different sizes");
    }

    std::vector<double> z(x.size());
    for (unsigned int i = 0; i < x.size(); i++) {
        z[i] = x[i] + scaling_factor * y[i];
    }
    return z;
}

double norm(std::vector<double>& x, std::string norm_value) {
    /* choose the right norm */
    if (norm_value == "inf") {
        return norm_inf(x);
    }
    else if (norm_value == "l2") {
        return norm_2(x);
    }
    else if (norm_value == "l2_squared") {
        return norm_2_squared(x);
    }
    else if (norm_value == "l1") {
        return norm_1(x);
    }
    else {
        throw std::out_of_range("The norm is not known");
    }
}

/* compute ||x||_1 */
double norm_1(std::vector<double>& x) {
    double norm = 0.;
    for (unsigned int i = 0; i < x.size(); i++) {
        norm += std::abs(x[i]);
    }
    return norm;
}

double norm_1(std::map<int, double>& x) {
    double norm = 0.;
    for (std::pair<int, double> term: x) {
        double xi = term.second;
        norm += std::abs(xi);
    }
    return norm;
}

// https://en.wikipedia.org/wiki/Matrix_norm#Special_cases
double norm_1(std::vector<std::map<int, double> >& m) {
    double norm = 0.;
    for (unsigned int j = 0; j < m.size(); j++) {
        double column_norm = norm_1(m[j]);
        norm = std::max(norm, column_norm);
    }
    return norm;
}

/* compute ||x||^2_2 */
double norm_2_squared(std::vector<double>& x) {
    double norm_squared = 0.;
    for (unsigned int i = 0; i < x.size(); i++) {
        norm_squared += x[i] * x[i];
    }
    return norm_squared;
}

/* compute ||x||_2 */
double norm_2(std::vector<double>& x) {
    return std::sqrt(norm_2_squared(x));
}

/* compute ||x||_infty */
double norm_inf(std::vector<double>& x, unsigned int length) {
    double norm = 0.;
    for (unsigned int i = 0; i < std::min<unsigned int>(length, x.size()); i++) {
        norm = std::max(norm, std::abs(x[i]));
    }
    return norm;
}

// https://en.wikipedia.org/wiki/Matrix_norm#Special_cases
double norm_inf(std::vector<std::map<int, double> >& m) {
    // compute maximum row index
    int number_rows = 0;
    for (unsigned int j = 0; j < m.size(); j++) {
        number_rows = std::max(number_rows, 1 + m[j].rbegin()->first);
    }
    // read the matrix column-wise and fill in the row_vectors norm vector
    std::vector<double> row_vectors(number_rows);
    for (unsigned int j = 0; j < m.size(); j++) {
        for (std::pair<const int, double> element: m[j]) {
            int i = element.first;
            double value = element.second;
            row_vectors[i] += std::abs(value);
        }
    }
    // compute the maximal component of the row_vectors vector
    double norm = 0.;
    for (unsigned int i = 0; i < row_vectors.size(); i++) {
        norm = std::max(norm, row_vectors[i]);
    }
    return norm;
}

double dot(std::vector<double>& x, std::vector<double>& y) {
    double dot = 0.;
    for (unsigned int i = 0; i < std::min(x.size(), y.size()); i++) {
        dot += x[i] * y[i];
    }
    return dot;
}

double dot(std::vector<double>& x, std::map<int, double>& y) {
    double dot = 0.;
    for (std::pair<int, double> term: y) {
        unsigned int i = term.first;
        double yi = term.second;
        if (i < x.size()) {
            dot += x[i] * yi;
        }
    }
    return dot;
}

double dot(std::map<int, double>& x, std::map<int, double>& y) {
    double dot = 0.;
    for (std::pair<int, double> term: x) {
        int i = term.first;
        double xi = term.second;
        try {
            dot += xi * y.at(i);
        }
        catch (std::out_of_range) {}
    }
    return dot;
}