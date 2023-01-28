// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAGRANGIANGRADIENT_H
#define UNO_LAGRANGIANGRADIENT_H

#include <iostream>
#include <vector>
#include "linear_algebra/Vector.hpp"

// Gradient of the Lagrangian
// Keep the objective and constraint contributions separate. This helps:
// - computing the KKT and FJ stationarity conditions
// - forming quasi-Newton matrices with an objective multiplier
template <typename T>
class LagrangianGradient {
public:
   std::vector<T> objective_contribution{};
   std::vector<T> constraints_contribution{};

   using value_type = T;

   explicit LagrangianGradient(size_t number_variables);
   [[nodiscard]] size_t size() const;
   [[nodiscard]] T operator[](size_t i) const;
   void resize(size_t new_number_variables);
};

template <typename T>
LagrangianGradient<T>::LagrangianGradient(size_t number_variables) :
      objective_contribution(number_variables),
      constraints_contribution(number_variables) {
}

template <typename T>
size_t LagrangianGradient<T>::size() const {
   return this->objective_contribution.size();
}

// access i-th element
template <typename T>
T LagrangianGradient<T>::operator[](size_t i) const {
   return this->objective_contribution[i] + this->constraints_contribution[i];
}

template <typename T>
void LagrangianGradient<T>::resize(size_t new_number_variables) {
   this->objective_contribution.resize(new_number_variables);
   this->constraints_contribution.resize(new_number_variables);
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const LagrangianGradient<T>& gradient) {
   for (size_t i: range(gradient.constraints_contribution.size())) {
      stream << gradient[i] << ' ';
   }
   stream << '\n';
   return stream;
}

#endif // UNO_LAGRANGIANGRADIENT_H


