// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOR_H
#define UNO_VECTOR_H

#include <iostream>
#include <string>
#include <vector>
#include <initializer_list>
#include "BLASVector.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   template <typename T>
   class Vector: public BLASVector<T> {
   public:
      using value_type = T;
      // iterators
      using iterator = typename std::vector<T>::iterator;
      using const_iterator = typename std::vector<T>::const_iterator;

      // constructors and destructor
      explicit Vector(size_t capacity = 0): vector(capacity) { }
      Vector(size_t capacity, T value): vector(capacity, value) { }
      Vector(std::initializer_list<T> initializer_list): vector(initializer_list) { }
      Vector(const Vector<T>& other): vector(other.vector) { }
      Vector(Vector<T>&& other) noexcept: vector(std::move(other.vector)) { }
      ~Vector() override = default;

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

      // operators
      using BLASVector<value_type>::operator=;
      using BLASVector<value_type>::operator+=;
      using BLASVector<value_type>::operator-=;

      [[nodiscard]] size_t size() const override {
         return this->vector.size();
      }

      T* data() override {
         return this->vector.data();
      }

      const T* data() const override {
         return this->vector.data();
      }

      T& operator[](size_t index) override {
         return this->vector[index];
      }

      const T& operator[](size_t index) const override {
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
         for (size_t index: Range(this->size())) {
            this->vector[index] = value;
         }
      }

      void print(std::ostream& stream) const {
         for (const T& element: *this) {
            stream << element << ' ';
         }
      }

   protected:
      std::vector<T> vector;
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