#include "Constraint.hpp"

Multipliers::Multipliers(int number_variables, int number_constraints): lower_bounds(number_variables), upper_bounds(number_variables), constraints(number_constraints) {
}
