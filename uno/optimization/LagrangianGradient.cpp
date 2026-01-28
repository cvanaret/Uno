// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <ostream>
#include "LagrangianGradient.hpp"

namespace uno {
   LagrangianGradient::LagrangianGradient(size_t number_variables) :
         objective_contribution(number_variables),
         constraints_contribution(number_variables) {
   }
   
   size_t LagrangianGradient::size() const {
      return this->objective_contribution.size();
   }

   // access i-th element
   double LagrangianGradient::operator[](size_t variable_index) const {
      return this->objective_contribution[variable_index] + this->constraints_contribution[variable_index];
   }

   void LagrangianGradient::resize(size_t number_variables) {
      this->objective_contribution.resize(number_variables);
      this->constraints_contribution.resize(number_variables);
   }
   
   std::ostream& operator<<(std::ostream& stream, const LagrangianGradient& gradient) {
      for (size_t variable_index: Range(gradient.constraints_contribution.size())) {
         stream << gradient[variable_index] << ' ';
      }
      stream << '\n';
      return stream;
   }
} // namespace