// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LagrangianGradient.hpp"

LagrangianGradient::LagrangianGradient(size_t number_variables) :
      objective_contribution(number_variables),
      constraints_contribution(number_variables) {
}

size_t LagrangianGradient::size() const {
   return this->objective_contribution.size();
}

// access i-th element
double LagrangianGradient::operator[](size_t i) const {
   return this->objective_contribution[i] + this->constraints_contribution[i];
}

/*
// compute ||x||_1
double LagrangianGradient::norm_1() const {
   double norm = 0.;
   for (size_t i: Range(this->objective_contribution.size())) {
      norm += std::abs(this->operator[](i));
   }
   return norm;
}
*/

void LagrangianGradient::resize(size_t new_number_variables) {
   this->objective_contribution.resize(new_number_variables);
   this->constraints_contribution.resize(new_number_variables);
}

std::ostream& operator<<(std::ostream& stream, const LagrangianGradient& gradient) {
   for (size_t i: Range(gradient.constraints_contribution.size())) {
      stream << gradient[i] << ' ';
   }
   stream << '\n';
   return stream;
}