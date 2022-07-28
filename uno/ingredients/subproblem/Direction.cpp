// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Direction.hpp"
#include "tools/Logger.hpp"
#include "linear_algebra/Vector.hpp"

Direction::Direction(size_t max_number_variables, size_t max_number_constraints):
   number_variables(max_number_variables), number_constraints(max_number_constraints),
   primals(max_number_variables), multipliers(max_number_variables, max_number_constraints) {
}

void Direction::set_dimensions(size_t number_variables, size_t number_constraints) {
   this->number_variables = number_variables;
   this->number_constraints = number_constraints;
}

std::string status_to_string(Status status) {
   switch (status) {
      case Status::OPTIMAL:
         return "optimal";
      case Status::UNBOUNDED_PROBLEM:
         return "unbounded problem";
      case Status::INFEASIBLE:
         return "infeasible";
      default:
         return "unknown status, something went wrong";
   }
}

std::ostream& operator<<(std::ostream& stream, const Direction& direction) {
   stream << "\nDirection:\n";
   stream << "d^* = ";
   print_vector(stream, direction.primals, 0, direction.number_variables);

   stream << "Status: " << status_to_string(direction.status) << '\n';

   stream << "objective = " << direction.objective << '\n';
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

   stream << "objective multiplier = " << direction.objective_multiplier << '\n';
   stream << "lower bound multipliers = ";
   print_vector(stream, direction.multipliers.lower_bounds);
   stream << "upper bound multipliers = ";
   print_vector(stream, direction.multipliers.upper_bounds);
   stream << "constraint multipliers = ";
   print_vector(stream, direction.multipliers.constraints);

   return stream;
}

ConstraintPartition::ConstraintPartition(size_t number_constraints) {
   this->feasible.reserve(number_constraints);
   this->infeasible.reserve(number_constraints);
}