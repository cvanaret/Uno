#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>
#include "Logger.hpp"
#include "SparseGradient.hpp"

std::vector<double> add_vectors(const std::vector<double>& x, const std::vector<double>& y, double scaling_factor = 1.);

double norm_1(const std::vector<double>& x);
double norm_1(const SparseGradient& x);
double norm_1(const std::vector<SparseGradient>& m);

double norm_2_squared(const std::vector<double>& x);
double norm_2_squared(const SparseGradient& x);
double norm_2(const std::vector<double>& x);
double norm_2(const SparseGradient& x);

double norm_inf(const std::vector<double>& x, unsigned int length = std::numeric_limits<unsigned int>::max());
double norm_inf(const SparseGradient& x);
double norm_inf(const std::vector<SparseGradient>& m);

double dot(const std::vector<double>& x, const std::vector<double>& y);
double dot(const std::vector<double>& x, const SparseGradient& y);
double dot(const SparseGradient& x, const SparseGradient& y);

template <typename T>
double norm(const T& x, const std::string& norm_value) {
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

template <typename T>
void print_vector(std::ostream &stream, const std::vector<T>& x, const char end='\n', unsigned int start = 0, unsigned int length = std::numeric_limits<unsigned int>::max()) {
    for (size_t i = start; i < std::min<unsigned int>(start + length, x.size()); i++) {
        stream << x[i] << " ";
    }
    stream << end;
    return;
}

template <typename T>
void print_vector(const Level& level, const std::vector<T>& x, const char end='\n', unsigned int start = 0, unsigned int length = std::numeric_limits<unsigned int>::max()) {
    for (size_t i = start; i < std::min<unsigned int>(start + length, x.size()); i++) {
        level << x[i] << " ";
    }
    level << end;
    return;
}

template <typename T>
void print_vector(const Level& level, const std::set<T>& x, const char end='\n') {
    for (T xi: x) {
        level << xi << " ";
    }
    level << end;
    return;
}

template <typename T, typename U>
void print_vector(std::ostream &stream, const std::map<T, U>& x, const char end='\n') {
    for (std::pair<T, U> element : x) {
        T i = element.first;
        U xi = element.second;
        stream << "x[" << i << "] = " << xi << ", ";
    }
    stream << end;
    return;
}

template <typename T, typename U>
void print_vector(const Level& level, const std::map<T, U>& x, const char end='\n') {
    for (std::pair<T, U> element : x) {
        T i = element.first;
        U xi = element.second;
        level << "x[" << i << "] = " << xi << ", ";
    }
    level << end;
    return;
}

template <typename T, typename U>
void print_vector(std::ostream &stream, const std::unordered_map<T, U>& x, const char end='\n') {
    for (std::pair<T, U> element : x) {
        T i = element.first;
        U xi = element.second;
        stream << "x[" << i << "] = " << xi << ", ";
    }
    stream << end;
    return;
}

template <typename T, typename U>
void print_vector(const Level& level, const std::unordered_map<T, U>& x, const char end='\n') {
    for (std::pair<T, U> element : x) {
        T i = element.first;
        U xi = element.second;
        level << "x[" << i << "] = " << xi << ", ";
    }
    level << end;
    return;
}

std::string join(std::vector<std::string> vector, const char* separator);

#endif // UTILS_H
