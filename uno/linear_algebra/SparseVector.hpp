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
   this->indices.push_back(index);
   this->values.push_back(value);
   this->number_nonzeros++;
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
double dot(const SparseVector& x, const SparseVector& y);
void print_vector(std::ostream &stream, const SparseVector& x, const char end);
void print_vector(const Level& level, const SparseVector& x, const char end);

/*
class SparseVector {
public:
   SparseVector(int capacity) {
      this->indices.reserve(capacity);
      this->elements.reserve(capacity);
      //std::cout << " *** Initialized with capacity " << capacity << "\n";
   }

   SparseVector() {
   }

   void reserve(size_t capacity) {
      this->indices.reserve(capacity);
      this->elements.reserve(capacity);
   }

   size_t size() const {
      return this->number_nonzeros;
   }

   double& operator[](size_t index) {
      size_t i = 0;
      while (i < this->number_nonzeros) {
         if (this->indices[i] == index) {
            return this->elements[i];
         }
         i++;
      }
      assert(this->number_nonzeros <= this->indices.capacity());
      this->indices[this->number_nonzeros] = index;
      double& result = this->elements[this->number_nonzeros];
      this->number_nonzeros++;
      return result;
   }

   const double& operator[](size_t index) const {
      size_t i = 0;
      while (i < this->number_nonzeros) {
         if (this->indices[i] == index) {
            return this->elements[i];
         }
         i++;
      }
      assert(false && "Entry not found");
   }

   double at(size_t i) const {
      return (*this)[i];
   }

   void erase(size_t index) {
      size_t i = 0;
      while (i < this->number_nonzeros) {
         if (this->indices[i] == index) {
            const int nnz = this->number_nonzeros;
            // swap with the last elements
            this->indices[i] = this->indices[nnz - 1];
            this->elements[i] = this->elements[nnz - 1];
            this->number_nonzeros--;
            return;
         }
         i++;
      }
   }

   void clear() {
      this->number_nonzeros = 0;
   }

   class iterator {
   public:
      iterator(const SparseVector& v, int i) : vector(&v) {
         current_index = i;
      }
      bool operator!=(const iterator x) const {
         return current_index != x.current_index;
      };
      iterator& operator++() {
         ++current_index;
         return *this;
      };
      std::pair<size_t, double> operator*() const {
         return std::pair(this->vector->indices[this->current_index], this->vector->elements[this->current_index]);
      };

   private:
      const SparseVector* vector;
      int current_index;
   };

   iterator begin() {
      return iterator(*this, 0);
   };
   iterator end() {
      return iterator(*this, this->number_nonzeros);
   };
   iterator begin() const {
      return iterator(*this, 0);
   };
   iterator end() const {
      return iterator(*this, this->number_nonzeros);
   };

//protected:
   std::vector<size_t> indices;
   std::vector<double> elements;
   size_t number_nonzeros{0};
};
*/

#endif // SPARSEVECTOR_H
