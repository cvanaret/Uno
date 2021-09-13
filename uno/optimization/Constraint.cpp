#include "Constraint.hpp"

ConstraintPartition::ConstraintPartition(size_t number_constraints) : constraint_feasibility(number_constraints) {
   this->feasible.reserve(number_constraints);
   this->infeasible.reserve(number_constraints);
}

Multipliers::Multipliers(size_t number_variables, size_t number_constraints) : lower_bounds(number_variables),
      upper_bounds(number_variables), constraints(number_constraints) {
}