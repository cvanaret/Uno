// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <ostream>
#include "Iterate.hpp"
#include "tools/Symbols.hpp"

namespace uno {
   Iterate::Iterate(size_t number_variables, size_t number_constraints) :
         number_variables(number_variables), number_constraints(number_constraints),
         primals(number_variables), multipliers(number_variables, number_constraints),
         residuals(number_variables) {
   }

   void Iterate::set_number_variables(size_t new_number_variables) {
      this->number_variables = new_number_variables;
      this->primals.resize(new_number_variables);
      this->residuals.lagrangian_gradient.resize(new_number_variables);
   }

   std::ostream& operator<<(std::ostream& stream, const Iterate& iterate) {
      stream << "Primal variables: " << iterate.primals << '\n';
      stream << "            " << symbols::top_pipe << " Constraint: " << iterate.multipliers.constraints << '\n';
      stream << "Multipliers " << symbols::pipe << " Lower bound: " << iterate.multipliers.lower_bounds << '\n';
      stream << "            " << symbols::bottom_pipe << " Upper bound: " << iterate.multipliers.upper_bounds << '\n';
      stream << "Primal feasibility: " << iterate.primal_feasibility << '\n';

      stream << "          " << symbols::top_pipe << " Stationarity: " << iterate.residuals.stationarity << '\n';
      stream << "Residuals " << symbols::pipe << " Complementarity: " << iterate.residuals.complementarity << '\n';
      stream << "          " << symbols::bottom_pipe << " Lagrangian gradient: " << iterate.residuals.lagrangian_gradient << '\n';

      stream << "                  " << symbols::top_pipe << " Infeasibility: " << iterate.progress.infeasibility << '\n';
      stream << "Progress measures " << symbols::pipe << " Optimality: " << iterate.progress.objective(1.) << '\n';
      stream << "                  " << symbols::bottom_pipe << " Auxiliary terms: " << iterate.progress.auxiliary << '\n';

      return stream;
   }
} // namespace