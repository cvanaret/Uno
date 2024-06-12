// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Direction.hpp"
#include "tools/Logger.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/VectorView.hpp"

Direction::Direction(size_t number_variables, size_t number_constraints) :
      number_variables(number_variables), number_constraints(number_constraints),
      primals(number_variables), multipliers(number_variables, number_constraints), feasibility_multipliers(number_variables, number_constraints),
      active_bounds(number_variables) {
}

void Direction::set_dimensions(size_t new_number_variables, size_t new_number_constraints) {
   this->number_variables = new_number_variables;
   this->number_constraints = new_number_constraints;
}

void Direction::reset() {
   this->primals.fill(0.);
   this->multipliers.reset();
   this->feasibility_multipliers.reset();
   this->active_bounds.at_lower_bound.clear();
   this->active_bounds.at_upper_bound.clear();
}

std::string status_to_string(SubproblemStatus status) {
   switch (status) {
      case SubproblemStatus::OPTIMAL:
         return "optimal subproblem";
      case SubproblemStatus::UNBOUNDED_PROBLEM:
         return "unbounded subproblem";
      case SubproblemStatus::INFEASIBLE:
         return "infeasible subproblem";
      default:
         return "unknown status, something went wrong";
   }
}

std::ostream& operator<<(std::ostream& stream, const Direction& direction) {
   stream << "Direction:\n";
   stream << "│ status: " << status_to_string(direction.status) << '\n';

   stream << "│ primals = "; print_vector(stream, view(direction.primals, 0, direction.number_variables));
   stream << "│ constraint multipliers = "; print_vector(stream, direction.multipliers.constraints);
   stream << "│ lower bound multipliers = "; print_vector(stream, direction.multipliers.lower_bounds);
   stream << "│ upper bound multipliers = "; print_vector(stream, direction.multipliers.upper_bounds);
   stream << "│ feasibility constraint multipliers = "; print_vector(stream, direction.feasibility_multipliers.constraints);
   stream << "│ feasibility lower bound multipliers = "; print_vector(stream, direction.feasibility_multipliers.lower_bounds);
   stream << "│ feasibility upper bound multipliers = "; print_vector(stream, direction.feasibility_multipliers.upper_bounds);

   stream << "│ objective = " << direction.subproblem_objective << '\n';
   stream << "│ norm = " << direction.norm << '\n';

   stream << "│ bound constraints active at lower bound =";
   for (size_t variable_index: direction.active_bounds.at_lower_bound) {
      stream << " x" << variable_index;
   }
   stream << '\n';
   stream << "│ bound constraints active at upper bound =";
   for (size_t variable_index: direction.active_bounds.at_upper_bound) {
      stream << " x" << variable_index;
   }
   stream << '\n';
   return stream;
}