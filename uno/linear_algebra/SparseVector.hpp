#ifndef SPARSEVECTOR_H
#define SPARSEVECTOR_H

#include <cassert>
#include <unordered_map>
#include <functional>
#include <ostream>
#include "tools/Logger.hpp"

using SparseVector = std::unordered_map<size_t, double>;

template <typename T>
class SparseVector2 {
public:
   SparseVector2(size_t capacity);
   void for_each(const std::function<void (size_t, T)>& f) const;
   size_t size() const;
   void reserve(size_t capacity);

   void insert(size_t index, T value);
   void transform(const std::function<T (T)>& f);
   void clear();

   template <typename U>
   friend std::ostream& operator<<(std::ostream& stream, const SparseVector2<U>& x);

protected:
   std::vector<size_t> indices;
   std::vector<T> values;
   size_t number_nonzeros{0};
};

// SparseVector2 methods
template <typename T>
SparseVector2<T>::SparseVector2(size_t capacity) {
   this->indices.reserve(capacity);
   this->values.reserve(capacity);
}

template <typename T>
void SparseVector2<T>::for_each(const std::function<void (size_t, T)>& f) const {
   for (size_t i = 0; i < this->number_nonzeros; i++) {
      f(this->indices[i], this->values[i]);
   }
}

template <typename T>
size_t SparseVector2<T>::size() const {
   return this->number_nonzeros;
}

template <typename T>
void SparseVector2<T>::reserve(size_t capacity) {
   this->indices.reserve(capacity);
   this->values.reserve(capacity);
}

template <typename T>
void SparseVector2<T>::insert(size_t index, T value) {
   auto element = std::find(begin(this->indices), end(this->indices), index);
   // if index is not found, add the new term, otherwise update it
   if (element == end(this->indices)) {
      this->indices.push_back(index);
      this->values.push_back(value);
      this->number_nonzeros++;
   }
   else {
      // *element is the index at which the index was found
      this->values[*element] += value;
   }
}

template <typename T>
void SparseVector2<T>::clear() {
   this->indices.clear();
   this->values.clear();
   this->number_nonzeros = 0;
}

template <typename T>
void SparseVector2<T>::transform(const std::function<T (T)>& f) {
   for (size_t i = 0; i < this->number_nonzeros; i++) {
      this->values[i] = f(this->values[i]);
   }
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const SparseVector2<T>& x) {
   stream << x.size() << " non zeros\n";
   x.for_each([&](size_t index, T entry) {
      stream << "index: " << index << " = " << entry << "\n";
   });
   return stream;
}

// free functions

void clear(SparseVector& x);
void scale(SparseVector& x, double scaling_factor);
double norm_1(const SparseVector2<double>& x);
// https://en.wikipedia.org/wiki/Matrix_norm#Special_cases
//double norm_1(const std::vector<SparseVector>& m);
double norm_2_squared(const SparseVector& x);
double norm_2(const SparseVector& x);
double norm_inf(const SparseVector& x);
// https://en.wikipedia.org/wiki/Matrix_norm#Special_cases
double norm_inf(const std::vector<SparseVector>& m);
double dot(const std::vector<double>& x, const SparseVector& y);
double dot(const std::vector<double>& x, const SparseVector2<double>& y);
double dot(const SparseVector& x, const SparseVector& y);
void print_vector(std::ostream &stream, const SparseVector& x, const char end);
void print_vector(const Level& level, const SparseVector& x, const char end);
void scale(SparseVector2<double>& x, double factor);

#endif // SPARSEVECTOR_H
