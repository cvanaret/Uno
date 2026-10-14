// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NORM_H
#define UNO_NORM_H

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include "symbolic/Range.hpp"

namespace uno {
   // norms of any array with elements of any type, or of dematerialized vectors given as (indices, component function)

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

   //**********
   // traits //
   //**********
   // materialized array: has size() and operator[](size_t)
   template <typename T, typename = void>
   struct is_indexable: std::false_type {};

   template <typename T>
   struct is_indexable<T, std::void_t<decltype(std::declval<const T&>()[std::declval<size_t>()]),
      decltype(std::declval<const T&>().size())>>: std::true_type {};

   template <typename... Ts>
   inline constexpr bool are_indexable_v = (is_indexable<Ts>::value && ...);

   // component function: i -> x[i]
   template <typename F>
   inline constexpr bool is_component_function_v = std::is_invocable_r_v<double, const F&, size_t>;

   //***********
   // l1 norm //
   //***********
   template <typename Array>
   std::enable_if_t<is_indexable<Array>::value, double> norm_1(const Array& x) {
      double result = 0.;
      for (size_t index: Range(x.size())) {
         result += std::abs(x[index]);
      }
      return result;
   }

   // l1 norm of several arrays
   template <typename Array, typename... Arrays>
   std::enable_if_t<are_indexable_v<Array, Arrays...>, double> norm_1(const Array& x, const Arrays& ... other_arrays) {
      return norm_1(x) + norm_1(other_arrays...);
   }

   // l1 norm of several dematerialized vectors: (indices, component_function)...
   template <typename Indices, typename ComponentFunction, typename... Rest>
   std::enable_if_t<is_component_function_v<ComponentFunction>, double> norm_1(const Indices& indices,
         const ComponentFunction& component_function, const Rest& ... rest) {
      static_assert(sizeof...(Rest) % 2 == 0, "Expected (indices, component_function) pairs");
      double result = 0.;
      for (const size_t index: indices) {
         result += std::abs(component_function(index));
      }
      if constexpr (sizeof...(Rest) > 0) {
         result += norm_1(rest...);
      }
      return result;
   }

   //*******************
   // l2 squared norm //
   //*******************
   template <typename Array>
   std::enable_if_t<is_indexable<Array>::value, double> norm_2_squared(const Array& x) {
      double result = 0.;
      for (size_t index: Range(x.size())) {
         result += x[index] * x[index];
      }
      return result;
   }

   // l2 squared norm of several arrays
   template <typename Array, typename... Arrays>
   std::enable_if_t<are_indexable_v<Array, Arrays...>, double> norm_2_squared(const Array& x, const Arrays& ... other_arrays) {
      return norm_2_squared(x) + norm_2_squared(other_arrays...);
   }

   // l2 squared norm of several dematerialized vectors: (indices, component_function)...
   template <typename Indices, typename ComponentFunction, typename... Rest>
   std::enable_if_t<is_component_function_v<ComponentFunction>, double> norm_2_squared(const Indices& indices,
         const ComponentFunction& component_function, const Rest& ... rest) {
      static_assert(sizeof...(Rest) % 2 == 0, "Expected (indices, component_function) pairs");
      double result = 0.;
      for (const size_t index: indices) {
         const double element = component_function(index);
         result += element * element;
      }
      if constexpr (sizeof...(Rest) > 0) {
         result += norm_2_squared(rest...);
      }
      return result;
   }

   //***********
   // l2 norm //
   //***********
   // arrays or (indices, component_function) pairs: dispatched by norm_2_squared
   template <typename... Arguments>
   double norm_2(const Arguments& ... arguments) {
      return std::sqrt(norm_2_squared(arguments...));
   }

   //**************
   // l_inf norm //
   //**************
   // NaN-propagating max
   inline double nan_max(double a, double b) {
      return (std::isnan(b) || a < b) ? b : a;
   }

   template <typename Array>
   std::enable_if_t<is_indexable<Array>::value, double> norm_inf(const Array& x) {
      double result = 0.;
      for (size_t index: Range(x.size())) {
         result = nan_max(result, std::abs(x[index]));
         // NaN is sticky: abs > NaN is always false afterwards
      }
      return result;
   }

   // inf norm of several arrays
   template <typename Array, typename... Arrays>
   std::enable_if_t<are_indexable_v<Array, Arrays...>, double> norm_inf(const Array& x, const Arrays& ... other_arrays) {
      return nan_max(norm_inf(x), norm_inf(other_arrays...));
   }

   // inf norm of several dematerialized vectors: (indices, component_function)...
   template <typename Indices, typename ComponentFunction, typename... Rest>
   std::enable_if_t<is_component_function_v<ComponentFunction>, double> norm_inf(const Indices& indices,
         const ComponentFunction& component_function, const Rest& ... rest) {
      static_assert(sizeof...(Rest) % 2 == 0, "Expected (indices, component_function) pairs");
      double result = 0.;
      for (const size_t index: indices) {
         result = nan_max(result, std::abs(component_function(index)));
      }
      if constexpr (sizeof...(Rest) > 0) {
         result = nan_max(result, norm_inf(rest...));
      }
      return result;
   }

   //*********************
   // dispatch function //
   //*********************
   // norm of at least one array, or of at least one (indices, component_function) pair
   template <typename... Arguments>
   double norm(Norm norm, const Arguments& ... arguments) {
      static_assert(sizeof...(Arguments) > 0, "norm needs at least one argument");
      if (norm == Norm::L1) {
         return norm_1(arguments...);
      }
      else if (norm == Norm::L2) {
         return norm_2(arguments...);
      }
      else if (norm == Norm::L2_SQUARED) {
         return norm_2_squared(arguments...);
      }
      else if (norm == Norm::INF) {
         return norm_inf(arguments...);
      }
      throw std::invalid_argument("The norm is not known");
   }
} // namespace

#endif // UNO_NORM_H