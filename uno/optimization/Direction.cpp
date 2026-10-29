// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
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
      this->status = SubproblemStatus::OPTIMAL;
   }

   bool Direction::is_finite() const {
      const auto finite = [](const Vector<double>& vector) {
         return std::all_of(vector.begin(), vector.end(), [](double value) { return std::isfinite(value); });
      };
      return finite(this->primals) && finite(this->multipliers.constraints) && finite(this->multipliers.lower_bounds) &&
         finite(this->multipliers.upper_bounds);
   }

   std::ostream& operator<<(std::ostream& stream, const Direction& direction) {
      stream << "Direction:\n";
      stream << symbols::top_pipe << " status: " << Direction::status_to_string(direction.status) << '\n';
      stream << symbols::pipe << " primals = " << direction.primals << '\n';
      stream << symbols::pipe << " constraint multipliers = " << direction.multipliers.constraints << '\n';
      stream << symbols::pipe << " lower bound multipliers = " << direction.multipliers.lower_bounds << '\n';
      stream << symbols::pipe << " upper bound multipliers = " << direction.multipliers.upper_bounds << '\n';
      //stream << symbols::pipe << " objective = " << direction.subproblem_objective << '\n';
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
         case SubproblemStatus::ERROR:
            return "subproblem solver error";
         default:
            return "unknown status, something went wrong";
      }
   }
} // namespace