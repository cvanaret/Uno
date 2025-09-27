// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOR_H
#define UNO_VECTOR_H

#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <initializer_list>
#include "symbolic/Range.hpp"

namespace uno {
   template <typename ElementType>
   class Vector {
   public:
      using value_type = ElementType;
      // iterators
      using iterator = typename std::vector<ElementType>::iterator;
      using const_iterator = typename std::vector<ElementType>::const_iterator;

      // constructors and destructor
      explicit Vector(size_t capacity = 0): vector(capacity) { }
      Vector(size_t capacity, ElementType value): vector(capacity, value) { }
      Vector(std::initializer_list<ElementType> initializer_list): vector(initializer_list) { }
      Vector(const Vector<ElementType>& other) noexcept : vector(other.vector) { }
      Vector(Vector<ElementType>&& other) noexcept : vector(std::move(other.vector)) { }
      ~Vector() = default;

      // copy assignment operator
      Vector<ElementType>& operator=(const Vector<ElementType>& other) {
         if (&other != this) {
            this->vector = other.vector;
         }
         return *this;
      }

      // move assignment operator
      Vector<ElementType>& operator=(Vector<ElementType>&& other) noexcept {
         if (&other != this) {
            this->vector = std::move(other.vector);
         }
         return *this;
      }

      // assignment operator from some expression
      template <typename Expression>
      Vector<ElementType>& operator=(const Expression& expression) {
         static_assert(std::is_same_v<typename Expression::value_type, ElementType>);
         assert(expression.size() <= this->size() && "The expression is larger than the current vector");
         for (size_t index: Range(expression.size())) {
            this->vector[index] = expression[index];
         }
         return *this;
      }

      // sum operator
      template <typename Expression>
      Vector<ElementType>& operator+=(const Expression& expression) {
         for (size_t index: Range(this->size())) {
            this->vector[index] += expression[index];
         }
         return *this;
      }

      // random access
      ElementType& operator[](size_t index) { return this->vector[index]; }
      const ElementType& operator[](size_t index) const { return this->vector[index]; }

      // size
      [[nodiscard]] size_t size() const { return this->vector.size(); }
      [[nodiscard]] bool empty() const { return (this->size() == 0); }

      // size and capacity
      void reserve(size_t new_capacity) { this->vector.reserve(new_capacity); }
      void resize(size_t new_size) { this->vector.resize(new_size); }

      // iterators
      iterator begin() noexcept { return this->vector.begin(); }
      iterator end() noexcept { return this->vector.end(); }
      const_iterator begin() const noexcept { return this->vector.cbegin(); }
      const_iterator end() const noexcept { return this->vector.cend(); }

      // insertion
      void push_back(ElementType element) { this->vector.push_back(element); }
      void emplace_back(ElementType element) { this->vector.emplace_back(element); }

      void fill(ElementType value) {
         for (size_t index: Range(this->size())) {
            this->vector[index] = value;
         }
      }

      void scale(ElementType factor) {
         for (size_t index: Range(this->size())) {
            this->vector[index] *= factor;
         }
      }

      void operator*=(ElementType factor) {
         this->scale(factor);
      }

      ElementType* data() { return this->vector.data(); }
      const ElementType* data() const { return this->vector.data(); }

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
      for (size_t index: Range(result.size())) {
         result[index] -= expression[index];
      }
   }

   template <typename Container>
   std::string join(const Container& vector, const std::string& separator) {
      std::string result{};
      size_t index = 0;
      for (const auto& element: vector) {
         if (0 < index) {
            result += separator;
         }
         result += element;
         ++index;
      }
      return result;
   }
} // namespace

#endif // UNO_VECTOR_H