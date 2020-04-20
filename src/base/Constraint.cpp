#include "Constraint.hpp"

ConstraintPartition::ConstraintPartition(int number_constraints): constraint_feasibility(number_constraints) {
}

Multipliers::Multipliers(int number_variables, int number_constraints): lower_bounds(number_variables), upper_bounds(number_variables), constraints(number_constraints) {
}