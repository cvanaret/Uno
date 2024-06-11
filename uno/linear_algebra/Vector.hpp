// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOR_H
#define UNO_VECTOR_H

#include <iostream>
#include <limits>
#include <vector>
#include <initializer_list>
#include "symbolic/Range.hpp"

template <typename ElementType>
class Vector {
public:
   using value_type = ElementType;
   // iterators
   using iterator = typename std::vector<ElementType>::iterator;
   using const_iterator = typename std::vector<ElementType>::const_iterator;

   // constructors and destructor
   explicit Vector(size_t capacity = 0): vector(capacity) { }
   explicit Vector(size_t capacity, ElementType value): vector(capacity, value) { }
   Vector(std::initializer_list<ElementType> initializer_list): vector(initializer_list) { }
   Vector(const Vector& other) noexcept : vector(other.vector) { }
   Vector(Vector&& other) noexcept : vector(std::move(other.vector)) { }
   ~Vector() = default;

   // copy assignment operator
   Vector& operator=(const Vector& other) {
      for (size_t index = 0; index < this->size(); index++) {
         this->vector[index] = other[index];
      }
      return *this;
   }

   // assignment operator from an expression
   template <typename Expression>
   Vector& operator=(const Expression& expression) {
      static_assert(std::is_same_v<typename Expression::value_type, ElementType>);
      for (size_t index = 0; index < this->size(); index++) {
         this->vector[index] = expression[index];
      }
      return *this;
   }

   // move assignment operator
   Vector& operator=(Vector&& other) noexcept {
      if (&other != this) {
         this->vector = std::move(other.vector);
      }
      return *this;
   }

   // random access
   ElementType& operator[](size_t index) { return this->vector[index]; }
   const ElementType& operator[](size_t index) const { return this->vector[index]; }

   // size
   [[nodiscard]] size_t size() const { return this->vector.size(); }
   [[nodiscard]] bool empty() const { return (this->size() == 0); }
   void resize(size_t new_size) { this->vector.resize(new_size); }

   // iterators
   iterator begin() noexcept { return this->vector.begin(); }
   iterator end() noexcept { return this->vector.end(); }
   const_iterator begin() const noexcept { return this->vector.cbegin(); }
   const_iterator end() const noexcept { return this->vector.cend(); }

   void fill(ElementType value) {
      for (size_t index = 0; index < this->size(); index++) {
         this->vector[index] = value;
      }
   }

   void scale(ElementType factor) {
      for (size_t index = 0; index < this->size(); index++) {
         this->vector[index] *= factor;
      }
   }

   ElementType* data() { return this->vector.data(); }
   const ElementType* data() const { return this->vector.data(); }

   // sum operator
   template <typename Expression>
   Vector& operator+=(const Expression& expression) {
      for (size_t index = 0; index < this->size(); index++) {
         this->vector[index] += expression[index];
      }
      return *this;
   }

   void print(std::ostream& stream) const {
      for (const ElementType& element: *this) {
         stream << element << ' ';
      }
   }

protected:
   std::vector<ElementType> vector;
};

// use && to allow temporaries (such as std::cout or logger DEBUG, WARNING, etc)
template <typename Array, typename Stream>
void print_vector(Stream&& stream, const Array& x) {
   for (size_t index: Range(x.size())) {
      stream << x[index] << " ";
   }
   stream << '\n';
}

template <typename ElementType>
std::ostream& operator<<(std::ostream& stream, const Vector<ElementType>& vector) {
   vector.print(stream);
   return stream;
}

// subtract operator
template <typename ResultExpression, typename Expression>
void operator-=(ResultExpression&& result, const Expression& expression) {
   for (size_t index = 0; index < result.size(); index++) {
      result[index] -= expression[index];
   }
}

#endif // UNO_VECTOR_H
