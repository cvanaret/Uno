#include "Constraint.hpp"

ConstraintPartition::ConstraintPartition(size_t number_constraints) : constraint_feasibility(number_constraints) {
}

Multipliers::Multipliers(size_t number_variables, size_t number_constraints) : lower_bounds(number_variables),
      upper_bounds(number_variables), constraints(number_constraints), objective(1.) {
}

void Multipliers::copy_from(const Multipliers& multipliers) {
   // simply copy all multipliers
   for (size_t j = 0; j < this->constraints.size(); j++) {
      this->constraints[j] = multipliers.constraints[j];
   }
   for (size_t i = 0; i < this->lower_bounds.size(); i++) {
      this->lower_bounds[i] = multipliers.lower_bounds[i];
      this->upper_bounds[i] = multipliers.upper_bounds[i];
   }
   this->objective = multipliers.objective;
}