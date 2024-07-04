// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NORM_H
#define UNO_NORM_H

#include <cmath>
#include "symbolic/Range.hpp"

// norms of any array with elements of any type

enum class Norm {L1, L2, L2_SQUARED, INF};

inline Norm norm_from_string(const std::string& norm_string) {
   if (norm_string == "L1") {
      return Norm::L1;
   }
   else if (norm_string == "L2") {
      return Norm::L2;
   }
   else if (norm_string == "L2_squared") {
      return Norm::L2_SQUARED;
   }
   else if (norm_string == "INF") {
      return Norm::INF;
   }
   throw std::invalid_argument("The norm " + norm_string + " is not known");
}

// generic norm function for iterators that return [key, value] pairs
// https://stackoverflow.com/questions/38701475/how-to-overload-function-for-different-iterator-value-types-in-c
template <typename KeyValueIterable, typename AccumulationFunction, typename ElementType = typename KeyValueIterable::value_type,
      typename std::enable_if_t<not std::is_member_function_pointer<decltype(&KeyValueIterable::operator[])>::value, int> = 0>
ElementType generic_norm(const KeyValueIterable& x, const AccumulationFunction& accumulation_function) {
   ElementType result{0};
   for (const auto [_, element]: x) {
      accumulation_function(result, element);
   }
   return result;
}

// generic norm function for iterators that return the elements
template <typename Array, typename AccumulationFunction, typename ElementType = typename Array::value_type>
ElementType generic_norm(const Array& x, const AccumulationFunction& accumulation_function) {
   ElementType result{0};
   for (size_t index: Range(x.size())) {
      accumulation_function(result, x[index]);
   }
   return result;
}

//***********
// l1 norm //
//***********
template <typename ElementType>
void norm_1_accumulation(ElementType& result, ElementType element) {
   result += std::abs(element);
}

template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_1(const Array& x) {
   return generic_norm(x, norm_1_accumulation<ElementType>);
}

// l1 norm of several arrays
template <typename Array, typename... Arrays, typename ElementType = typename Array::value_type>
ElementType norm_1(const Array& x, const Arrays& ... other_arrays) {
   return norm_1(x) + norm_1(other_arrays...);
}

//*******************
// l2 squared norm //
//*******************
template <typename ElementType>
void norm_2_squared_accumulation(ElementType& result, ElementType element) {
   result += element * element;
}

template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_2_squared(const Array& x) {
   return generic_norm(x, norm_2_squared_accumulation<ElementType>);
}

// l2 squared norm of several arrays
template <typename Array, typename... Arrays, typename ElementType = typename Array::value_type>
ElementType norm_2_squared(const Array& x, const Arrays& ... other_arrays) {
   return norm_2_squared(x) + norm_2_squared(other_arrays...);
}

//***********
// l2 norm //
//***********
template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_2(const Array& x) {
   return std::sqrt(norm_2_squared(x));
}

// l2 norm of several arrays
template <typename Array, typename... Arrays>
typename Array::value_type norm_2(const Array& x, const Arrays& ... other_arrays) {
   return std::sqrt(norm_2_squared(x) + norm_2_squared(other_arrays...));
}

//**************
// l_inf norm //
//**************
template <typename ElementType>
void norm_inf_accumulation(ElementType& result, ElementType element) {
   result = std::max(result, std::abs(element));
}

template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_inf(const Array& x) {
   return generic_norm(x, norm_inf_accumulation<ElementType>);
}

// inf norm of several arrays
template <typename Array, typename... Arrays, typename ElementType = typename Array::value_type>
ElementType norm_inf(const Array& x, const Arrays& ... other_arrays) {
   return std::max(norm_inf(x), norm_inf(other_arrays...));
}

//*********************
// dispatch function //
//*********************
// norm of at least one array
template <typename Array, typename... Arrays, typename ElementType = typename Array::value_type>
ElementType norm(Norm norm, const Array& x, const Arrays& ... other_arrays) {
   if (norm == Norm::L1) {
      return norm_1(x, other_arrays...);
   }
   else if (norm == Norm::L2) {
      return norm_2(x, other_arrays...);
   }
   else if (norm == Norm::L2_SQUARED) {
      return norm_2_squared(x, other_arrays...);
   }
   else if (norm == Norm::INF) {
      return norm_inf(x, other_arrays...);
   }
   throw std::invalid_argument("The norm is not known");
}

#endif // UNO_NORM_H