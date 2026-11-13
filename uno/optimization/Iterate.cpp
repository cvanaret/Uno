// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <ostream>
#include "Iterate.hpp"
#include "tools/Symbols.hpp"

namespace uno {
   Iterate::Iterate(size_t number_variables, size_t number_constraints) :
         number_variables(number_variables), number_constraints(number_constraints),
         primals(number_variables), multipliers(number_variables, number_constraints),
         scaled_residuals(number_variables) {
   }

   Iterate::Iterate(): Iterate(0, 0) {}

   void Iterate::set_dimensions(size_t new_number_variables, size_t new_number_constraints) {
      this->number_variables = new_number_variables;
      this->number_constraints = new_number_constraints;
      this->primals.resize(new_number_variables);
      this->multipliers.constraints.resize(new_number_constraints);
      this->multipliers.lower_bounds.resize(new_number_variables);
      this->multipliers.upper_bounds.resize(new_number_variables);
      this->scaled_residuals.lagrangian_gradient.resize(new_number_variables);
   }

   std::ostream& operator<<(std::ostream& stream, const Iterate& iterate) {
      stream << "Primal variables: " << iterate.primals << '\n';
      stream << "            " << symbols::top_pipe << " Constraint: " << iterate.multipliers.constraints << '\n';
      stream << "Multipliers " << symbols::pipe << " Lower bound: " << iterate.multipliers.lower_bounds << '\n';
      stream << "            " << symbols::bottom_pipe << " Upper bound: " << iterate.multipliers.upper_bounds << '\n';

      stream << "          " << symbols::top_pipe << " Primal infeasibility: " << iterate.scaled_residuals.primal_infeasibility << '\n';
      stream << "Residuals " << symbols::pipe << " Stationarity: " << iterate.scaled_residuals.stationarity << '\n';
      stream << "          " << symbols::bottom_pipe << " Complementarity: " << iterate.scaled_residuals.complementarity << '\n';

      stream << "Lagrangian gradient: " << iterate.scaled_residuals.lagrangian_gradient << '\n';

      stream << "                  " << symbols::top_pipe << " Infeasibility: " << iterate.progress.infeasibility << '\n';
      stream << "Progress measures " << symbols::pipe << " Optimality: ";
      if (iterate.progress.objective != nullptr) {
         stream << iterate.progress.objective(1.) << '\n';
      }
      else {
         stream << "unset" << '\n';
      }
      stream << "                  " << symbols::bottom_pipe << " Auxiliary terms: " << iterate.progress.auxiliary << '\n';

      return stream;
   }
} // namespace