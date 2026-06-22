// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Direction.hpp"
#include "tools/Logger.hpp"
#include "tools/Symbols.hpp"

namespace uno {
   Direction::Direction(size_t number_variables, size_t number_constraints) :
         number_variables(number_variables), number_constraints(number_constraints),
         primals(number_variables), multipliers(number_variables, number_constraints) {
   }

   void Direction::set_dimensions(size_t new_number_variables, size_t new_number_constraints) {
      this->primals.resize(new_number_variables);
      this->multipliers.constraints.resize(new_number_constraints);
      this->multipliers.lower_bounds.resize(new_number_variables);
      this->multipliers.upper_bounds.resize(new_number_variables);
      this->number_variables = new_number_variables;
      this->number_constraints = new_number_constraints;
   }

   void Direction::reset() {
      this->primals.fill(0.);
      this->multipliers.reset();
   }

   std::ostream& operator<<(std::ostream& stream, const Direction& direction) {
      stream << "Direction:\n";
      stream << symbols::top_pipe << " status: " << Direction::status_to_string(direction.status) << '\n';
      stream << symbols::pipe << " primals = "; print_vector(stream, direction.primals);
      stream << symbols::pipe << " constraint multipliers = "; print_vector(stream, direction.multipliers.constraints);
      stream << symbols::pipe << " lower bound multipliers = "; print_vector(stream, direction.multipliers.lower_bounds);
      stream << symbols::pipe << " upper bound multipliers = "; print_vector(stream, direction.multipliers.upper_bounds);
      stream << symbols::pipe << " objective = " << direction.subproblem_objective << '\n';
      stream << symbols::bottom_pipe << " norm = " << direction.norm << '\n';
      return stream;
   }

   std::string Direction::status_to_string(SubproblemStatus status) {
      switch (status) {
         case SubproblemStatus::OPTIMAL:
            return "optimal";
         case SubproblemStatus::UNBOUNDED_PROBLEM:
            return "unbounded subproblem";
         case SubproblemStatus::INFEASIBLE:
            return "infeasible subproblem";
         default:
            return "unknown status, something went wrong";
      }
   }
} // namespace