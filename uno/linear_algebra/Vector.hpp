// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOR_H
#define UNO_VECTOR_H

#include <iostream>
#include <vector>
#include <cassert>
#include "tools/Logger.hpp"
#include "symbolic/Range.hpp"

// result <- x + scaling_factor * y
template <typename ElementType>
void add_vectors(const std::vector<ElementType>& x, const std::vector<ElementType>& y, ElementType scaling_factor, std::vector<ElementType>& result) {
   assert(x.size() <= y.size() && "Vector.add_vectors: x is longer than y");
   assert(x.size() <= result.size() && "Vector.add_vectors: result is not long enough");

   for (size_t index: Range(x.size())) {
      result[index] = x[index] + scaling_factor * y[index];
   }
}

// result <- a*x + b*y
template <typename ElementType>
void add_vectors(ElementType a, const std::vector<ElementType>& x, ElementType b, const std::vector<ElementType>& y, std::vector<ElementType>& result) {
   assert(x.size() <= y.size() && "Vector.add_vectors: x is longer than y");
   assert(x.size() <= result.size() && "Vector.add_vectors: result is not long enough");

   for (size_t index: Range(x.size())) {
      result[index] = a*x[index] + b*y[index];
   }
}

template <typename ElementType>
void initialize_vector(std::vector<ElementType>& x, ElementType value) {
   for (ElementType& element: x) {
      element = ElementType(value);
   }
}

template <typename ElementType>
void scale(std::vector<ElementType>& x, ElementType scaling_factor) {
   for (ElementType& xi: x) {
      xi *= scaling_factor;
   }
}

template <typename ElementType>
ElementType dot_product(const std::vector<ElementType>& x, const std::vector<ElementType>& y) {
   assert(x.size() == y.size() && "The vectors do not have the same size.");

   ElementType result = ElementType(0);
   for (size_t index: Range(x.size())) {
      result += x[index]*y[index];
   }
   return result;
}

template <typename ElementType>
void copy_from(std::vector<ElementType>& destination, const std::vector<ElementType>& source, size_t length = std::numeric_limits<size_t>::max()) {
   length = std::min(length, std::min(source.size(), destination.size()));
   const auto source_start_position = std::cbegin(source);
   const auto source_end_position = source_start_position + length;
   const auto destination_position = std::begin(destination);
   std::copy(source_start_position, source_end_position, destination_position);
}

// use && to allow temporaries (such as std::cout or logger DEBUG, WARNING, etc)
template <typename Array, typename Stream>
void print_vector(Stream&& stream, const Array& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max()) {
   for (size_t index: Range(start, std::min(start + length, x.size()))) {
      stream << x[index] << " ";
   }
   stream << '\n';
}

// check that an array of integers is in increasing order (x[i] <= x[i+1])
template <typename Array>
bool in_increasing_order(const Array& array, size_t length) {
   size_t index = 0;
   while (index < length - 1) {
      if (array[index] > array[index + 1]) {
         return false;
      }
      index++;
   }
   return true;
}

/*
// see here: https://stackoverflow.com/questions/10173623/override-operators-for-an-existing-class
template <typename T>
void operator+=(std::vector<T>& vector, const T& value) {
   for (T& element: vector) {
      element += value;
   }
}
*/

#endif // UNO_VECTOR_H
