#include <cmath>
#include "Constraint.hpp"

ConstraintPartition::ConstraintPartition(size_t number_constraints) {
   this->feasible.reserve(number_constraints);
   this->infeasible.reserve(number_constraints);
}

Multipliers::Multipliers(size_t number_variables, size_t number_constraints) : lower_bounds(number_variables),
      upper_bounds(number_variables), constraints(number_constraints) {
}

bool is_finite_lower_bound(double value) {
   return -std::numeric_limits<double>::infinity() < value;
}

bool is_finite_upper_bound(double value) {
   return value < std::numeric_limits<double>::infinity();
}