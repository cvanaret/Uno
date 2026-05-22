// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOR_H
#define UNO_VECTOR_H

#include <string>
#include <vector>
#include <initializer_list>
#include "VectorView.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   template <typename T>
   class Vector {
   public:
      using value_type = typename std::vector<T>::value_type;
      // iterators
      using iterator = typename std::vector<T>::iterator;
      using const_iterator = typename std::vector<T>::const_iterator;

      // constructors and destructor
      explicit Vector(size_t capacity = 0): vector(capacity) { }
      Vector(size_t capacity, T value): vector(capacity, value) { }
      Vector(std::initializer_list<T> initializer_list): vector(initializer_list) { }
      Vector(const Vector<T>& other): vector(other.vector) { }
      Vector(Vector<T>&& other) noexcept: vector(std::move(other.vector)) { }
      ~Vector() = default;

      // copy assignment operator
      Vector<T>& operator=(const Vector<T>& other) {
         if (&other != this) {
            this->vector = other.vector;
         }
         return *this;
      }

      // move assignment operator
      Vector<T>& operator=(Vector<T>&& other) noexcept {
         if (&other != this) {
            this->vector = std::move(other.vector);
         }
         return *this;
      }

      [[nodiscard]] size_t size() const {
         return this->vector.size();
      }

      [[nodiscard]] bool empty() const {
         return (this->size() == 0);
      }

      T* data() {
         return this->vector.data();
      }

      const T* data() const {
         return this->vector.data();
      }

      T& operator[](size_t index) {
         return this->vector[index];
      }

      const T& operator[](size_t index) const {
         return this->vector[index];
      }

      // size and capacity
      void reserve(size_t new_capacity) {
         this->vector.reserve(new_capacity);
      }

      void resize(size_t new_size) {
         this->vector.resize(new_size);
      }

      // iterators
      iterator begin() noexcept { return this->vector.begin(); }
      iterator end() noexcept { return this->vector.end(); }
      const_iterator begin() const noexcept { return this->vector.cbegin(); }
      const_iterator end() const noexcept { return this->vector.cend(); }

      // insertion
      void push_back(T element) {
         this->vector.push_back(element);
      }

      template <typename... Args>
      void emplace_back(Args&&... args) {
         this->vector.emplace_back(std::forward<Args>(args)...);
      }

      void fill(T value) {
         this->view().fill(value);
      }

      void print(std::ostream& stream) const {
         for (const T& element: *this) {
            stream << element << ' ';
         }
      }

      // mathematical operators: delegate to underlying VectorView

      template <typename Expression>
      Vector<T>& operator=(Expression&& expression) {
         this->view() = std::forward<Expression>(expression);
         return *this;
      }

      template <typename Expression>
      Vector<T>& operator+=(Expression&& expression) {
         this->view() += std::forward<Expression>(expression);
         return *this;
      }

      template <typename Expression>
      Vector<T>& operator-=(Expression&& expression) {
         this->view() -= std::forward<Expression>(expression);
         return *this;
      }

      void scale(T factor) {
         this->view().scale(factor);
      }

   protected:
      std::vector<T> vector;

      VectorView<const T> view() const noexcept {
         return uno::view(this->data(), this->size());
      }

      VectorView<T> view() noexcept {
         return uno::view(this->data(), this->size());
      }
   };

   // use && to allow temporaries (such as std::cout or logger DEBUG, WARNING, etc)
   template <typename Array, typename Stream>
   void print_vector(Stream&& stream, const Array& x) {
      for (size_t index: Range(x.size())) {
         stream << x[index] << " ";
      }
      stream << '\n';
   }

   template <typename T>
   std::ostream& operator<<(std::ostream& stream, const Vector<T>& vector) {
      vector.print(stream);
      return stream;
   }

   template <typename Container>
   std::string join(const Container& container, const std::string& separator) {
      std::string result{};
      size_t index = 0;
      for (const auto& element: container) {
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