// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <iostream>
#include <cassert>
#include <utility>
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Range.hpp"
#include "tools/Infinity.hpp"

std::map<FunctionType, std::string> Model::type_to_string = {
      {LINEAR, "linear"},
      {QUADRATIC, "quadratic"},
      {NONLINEAR, "nonlinear"}
};

// abstract Problem class
Model::Model(std::string name, size_t number_variables, size_t number_constraints, FunctionType type) :
      name(std::move(name)), number_variables(number_variables), number_constraints(number_constraints), problem_type(type),
      equality_constraints(this->number_constraints),
      inequality_constraints(this->number_constraints),
      linear_constraints(this->number_constraints),
      slacks(this->number_constraints) {
}

void Model::determine_bounds_types(std::vector<Interval>& bounds, std::vector<BoundType>& status) {
   assert(bounds.size() == status.size());
   // build the "status" vector as a mapping (map/transform operation) of the "bounds" vector
   std::transform(begin(bounds), end(bounds), begin(status), [](const Interval& bounds_i) {
      if (bounds_i.lb == bounds_i.ub) {
         return EQUAL_BOUNDS;
      }
      else if (is_finite(bounds_i.lb) && is_finite(bounds_i.ub)) {
         return BOUNDED_BOTH_SIDES;
      }
      else if (is_finite(bounds_i.lb)) {
         return BOUNDED_LOWER;
      }
      else if (is_finite(bounds_i.ub)) {
         return BOUNDED_UPPER;
      }
      else {
         return UNBOUNDED;
      }
   });
}

void Model::determine_constraints() {
   size_t current_equality_constraint = 0;
   size_t current_inequality_constraint = 0;
   for (size_t j = 0; j < this->number_constraints; j++) {
      if (this->get_constraint_bound_type(j) == EQUAL_BOUNDS) {
         this->equality_constraints.insert(j, current_equality_constraint);
         current_equality_constraint++;
      }
      else {
         this->inequality_constraints.insert(j, current_inequality_constraint);
         current_inequality_constraint++;
      }
   }
}

void Model::project_point_onto_bounds(std::vector<double>& x) const {
   for (size_t i = 0; i < x.size(); i++) {
      if (x[i] < this->get_variable_lower_bound(i)) {
         x[i] = this->get_variable_lower_bound(i);
      }
      else if (this->get_variable_upper_bound(i) < x[i]) {
         x[i] = this->get_variable_upper_bound(i);
      }
   }
}

bool Model::is_constrained() const {
   return (0 < this->number_constraints);
}

double Model::compute_constraint_lower_bound_violation(double constraint, size_t j) const {
   const double lower_bound = this->get_constraint_lower_bound(j);
   return std::max(0., lower_bound - constraint);
}

double Model::compute_constraint_upper_bound_violation(double constraint, size_t j) const {
   const double upper_bound = this->get_constraint_upper_bound(j);
   return std::max(0., constraint - upper_bound);
}

double Model::compute_constraint_violation(double constraint, size_t j) const {
   const double lower_bound_violation = this->compute_constraint_lower_bound_violation(constraint, j);
   const double upper_bound_violation = this->compute_constraint_upper_bound_violation(constraint, j);
   return std::max(lower_bound_violation, upper_bound_violation);
}

// compute ||c_S|| for a given set of constraints
double Model::compute_constraint_violation(const std::vector<double>& constraints, const std::vector<size_t>& constraint_set,
      Norm residual_norm) const {
   auto residual_function = [&](size_t j) {
      return this->compute_constraint_violation(constraints[j], j);
   };
   return norm(residual_function, constraint_set, residual_norm);
}

// compute ||c||
double Model::compute_constraint_violation(const std::vector<double>& constraints, Norm residual_norm) const {
   // create a lambda to avoid allocating an std::vector
   auto residual_function = [&](size_t j) {
      return this->compute_constraint_violation(constraints[j], j);
   };
   return norm(residual_function, Range(constraints.size()), residual_norm);
}

// complementary slackness error
double Model::compute_complementarity_error(const std::vector<double>& x, const std::vector<double>& constraints,
      const std::vector<double>& constraint_multipliers, const std::vector<double>& lower_bounds_multipliers,
      const std::vector<double>& upper_bounds_multipliers) const {
   // TODO use the objective multiplier
   double error = 0.;
   // bound constraints
   for (size_t i = 0; i < this->number_variables; i++) {
      if (is_finite(this->get_variable_lower_bound(i))) {
         error += std::abs(lower_bounds_multipliers[i] * (x[i] - this->get_variable_lower_bound(i)));
      }
      if (is_finite(this->get_variable_upper_bound(i))) {
         error += std::abs(upper_bounds_multipliers[i] * (x[i] - this->get_variable_upper_bound(i)));
      }
   }
   // constraints
   for (size_t j = 0; j < this->number_constraints; j++) {
      const double multiplier_j = constraint_multipliers[j];
      const double lower_bound = this->get_constraint_lower_bound(j);
      const double upper_bound = this->get_constraint_upper_bound(j);
      if (constraints[j] < lower_bound) { // violated lower
         // the optimal multiplier is 1
         error += std::abs((1. - multiplier_j) * (lower_bound - constraints[j]));
      }
      else if (upper_bound < constraints[j]) { // violated upper
         // the optimal multiplier is -1
         error += std::abs((1. + multiplier_j) * (constraints[j] - upper_bound));
      }
      else if (is_finite(lower_bound) && 0. < multiplier_j) { // lower bound
         error += std::abs(multiplier_j * (constraints[j] - lower_bound));
      }
      else if (is_finite(upper_bound) && multiplier_j < 0.) { // upper bound
         error += std::abs(multiplier_j * (constraints[j] - upper_bound));
      }
   }
   return error;
}
