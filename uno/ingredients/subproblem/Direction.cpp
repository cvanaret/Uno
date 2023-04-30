// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Direction.hpp"
#include "tools/Logger.hpp"
#include "linear_algebra/Vector.hpp"

Direction::Direction(size_t max_number_variables, size_t max_number_constraints) :
      number_variables(max_number_variables), number_constraints(max_number_constraints),
      primals(max_number_variables), multipliers(max_number_variables, max_number_constraints),
      active_set(max_number_variables, max_number_constraints) {
}

void Direction::set_dimensions(size_t new_number_variables, size_t new_number_constraints) {
   this->number_variables = new_number_variables;
   this->number_constraints = new_number_constraints;
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
   stream << "\nDirection:\n";
   stream << "Status: " << status_to_string(direction.status) << '\n';

   stream << "d^* = ";
   print_vector(stream, direction.primals, 0, direction.number_variables);
   stream << "constraint multipliers = ";
   print_vector(stream, direction.multipliers.constraints);
   stream << "lower bound multipliers = ";
   print_vector(stream, direction.multipliers.lower_bounds);
   stream << "upper bound multipliers = ";
   print_vector(stream, direction.multipliers.upper_bounds);
   stream << "objective multiplier = " << direction.objective_multiplier << '\n';

   stream << "objective = " << direction.subproblem_objective << '\n';
   stream << "norm = " << direction.norm << '\n';

   stream << "bound constraints active at lower bound =";
   for (size_t i: direction.active_set.bounds.at_lower_bound) {
      stream << " x" << i;
   }
   stream << '\n';
   stream << "bound constraints active at upper bound =";
   for (size_t i: direction.active_set.bounds.at_upper_bound) {
      stream << " x" << i;
   }
   stream << '\n';

   stream << "constraints at lower bound =";
   for (size_t j: direction.active_set.constraints.at_lower_bound) {
      stream << " c" << j;
   }
   stream << '\n';
   stream << "constraints at upper bound =";
   for (size_t j: direction.active_set.constraints.at_upper_bound) {
      stream << " c" << j;
   }
   stream << '\n';

   if (direction.constraint_partition.has_value()) {
      const ConstraintPartition& constraint_partition = direction.constraint_partition.value();
      stream << "general feasible =";
      for (size_t j: constraint_partition.feasible) {
         stream << " c" << j;
      }
      stream << "\ngeneral lower infeasible =";
      for (size_t j: constraint_partition.lower_bound_infeasible) {
         stream << " c" << j;
      }
      stream << "\ngeneral upper infeasible =";
      for (size_t j: constraint_partition.upper_bound_infeasible) {
         stream << " c" << j;
      }
      stream << '\n';
   }
   return stream;
}

ActiveConstraints::ActiveConstraints(size_t capacity) {
   this->at_lower_bound.reserve(capacity);
   this->at_upper_bound.reserve(capacity);
}

ActiveSet::ActiveSet(size_t number_variables, size_t number_constraints): constraints(number_constraints), bounds(number_variables) {
}

ConstraintPartition::ConstraintPartition(size_t number_constraints) {
   this->feasible.reserve(number_constraints);
   this->infeasible.reserve(number_constraints);
   this->lower_bound_infeasible.reserve(number_constraints);
   this->upper_bound_infeasible.reserve(number_constraints);
}