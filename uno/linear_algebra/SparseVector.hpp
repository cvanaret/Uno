// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPARSEVECTOR_H
#define UNO_SPARSEVECTOR_H

#include <cassert>
#include <functional>
#include "tools/Logger.hpp"

// SparseVector is a sparse vector that uses contiguous memory. It contains:
// - a vector of indices of type size_t
// - a vector of values of type T
// the indices are unique but not sorted
template <typename T>
class SparseVector {
public:
   explicit SparseVector(size_t capacity = 0);
   void for_each(const std::function<void (size_t, T)>& f) const;
   void for_each_index(const std::function<void (size_t)>& f) const;
   void for_each_value(const std::function<void (T)>& f) const;
   [[nodiscard]] size_t size() const;
   void reserve(size_t capacity);

   void insert(size_t index, T value);
   void transform(const std::function<T (T)>& f);
   void clear();
   [[nodiscard]] bool empty() const;

   template <typename U>
   friend std::ostream& operator<<(std::ostream& stream, const SparseVector<U>& x);

protected:
   std::vector<size_t> indices{};
   std::vector<T> values{};
   size_t number_nonzeros{0};
};

// SparseVector methods
template <typename T>
SparseVector<T>::SparseVector(size_t capacity) {
   this->reserve(capacity);
}

template <typename T>
void SparseVector<T>::for_each(const std::function<void (size_t, T)>& f) const {
   for (size_t i = 0; i < this->number_nonzeros; i++) {
      f(this->indices[i], this->values[i]);
   }
}

template <typename T>
void SparseVector<T>::for_each_index(const std::function<void(size_t)>& f) const {
   for (size_t i = 0; i < this->number_nonzeros; i++) {
      f(this->indices[i]);
   }
}

template <typename T>
void SparseVector<T>::for_each_value(const std::function<void(T)>& f) const {
   for (size_t i = 0; i < this->number_nonzeros; i++) {
      f(this->values[i]);
   }
}

template <typename T>
size_t SparseVector<T>::size() const {
   return this->number_nonzeros;
}

template <typename T>
void SparseVector<T>::reserve(size_t capacity) {
   this->indices.reserve(capacity);
   this->values.reserve(capacity);
}

template <typename T>
void SparseVector<T>::insert(size_t index, T value) {
   this->indices.push_back(index);
   this->values.push_back(value);
   this->number_nonzeros++;
}

template <typename T>
void SparseVector<T>::clear() {
   this->indices.clear();
   this->values.clear();
   this->number_nonzeros = 0;
}

template <typename T>
bool SparseVector<T>::empty() const {
   return (this->number_nonzeros == 0);
}

template <typename T>
void SparseVector<T>::transform(const std::function<T (T)>& f) {
   for (size_t i = 0; i < this->number_nonzeros; i++) {
      this->values[i] = f(this->values[i]);
   }
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const SparseVector<T>& x) {
   stream << "sparse vector with " << x.size() << " non zeros\n";
   x.for_each([&](size_t index, T entry) {
      stream << "index " << index << ", value " << entry << '\n';
   });
   return stream;
}

// free functions

double norm_1(const SparseVector<double>& x);
double norm_inf(const SparseVector<double>& x);
double dot(const std::vector<double>& x, const SparseVector<double>& y);
double dot(const std::vector<double>& x, const SparseVector<double>& y, const std::function<bool (size_t i)>& predicate);
void scale(SparseVector<double>& x, double factor);

#endif // UNO_SPARSEVECTOR_H
