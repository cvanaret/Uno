#include <cmath>
#include <iostream>
#include <cassert>
#include <utility>
#include "Problem.hpp"
#include "linear_algebra/Vector.hpp"

std::map<FunctionType, std::string> Problem::type_to_string = {
      {LINEAR, "linear"},
      {QUADRATIC, "quadratic"},
      {NONLINEAR, "nonlinear"}
};

// abstract Problem class
Problem::Problem(std::string name, size_t number_variables, size_t number_constraints, FunctionType type) :
      name(std::move(name)), number_variables(number_variables), number_constraints(number_constraints), problem_type(type),
      // allocate all vectors
      variables_bounds(number_variables), variable_status(number_variables),
      constraint_bounds(number_constraints), constraint_type(number_constraints), constraint_status(number_constraints) {
}

std::vector<double> Problem::evaluate_constraints(const std::vector<double>& x) const {
   // allocate an empty vector
   std::vector<double> constraints(this->number_constraints);
   // evaluate the constraints in place
   this->evaluate_constraints(x, constraints);
   return constraints;
}

void Problem::determine_bounds_types(std::vector<Range>& bounds, std::vector<ConstraintType>& status) {
   assert(bounds.size() == status.size());
   // build the "status" vector as a mapping (map/transform operation) of the "bounds" vector
   std::transform(begin(bounds), end(bounds), begin(status), [](const Range& bounds_i) {
      if (bounds_i.lb == bounds_i.ub) {
         return EQUAL_BOUNDS;
      }
      else if (is_finite_lower_bound(bounds_i.lb) && is_finite_upper_bound(bounds_i.ub)) {
         return BOUNDED_BOTH_SIDES;
      }
      else if (is_finite_lower_bound(bounds_i.lb)) {
         return BOUNDED_LOWER;
      }
      else if (is_finite_upper_bound(bounds_i.ub)) {
         return BOUNDED_UPPER;
      }
      else {
         return UNBOUNDED;
      }
   });
}

void Problem::determine_constraints() {
   size_t current_equality_constraint = 0;
   size_t current_inequality_constraint = 0;
   for (size_t j = 0; j < this->number_constraints; j++) {
      if (this->constraint_status[j] == EQUAL_BOUNDS) {
         this->equality_constraints[j] = current_equality_constraint;
         current_equality_constraint++;
      }
      else {
         this->inequality_constraints[j] = current_inequality_constraint;
         current_inequality_constraint++;
      }
   }
}

void Problem::project_point_in_bounds(std::vector<double>& x) const {
   for (size_t i = 0; i < x.size(); i++) {
      if (x[i] < this->variables_bounds[i].lb) {
         x[i] = this->variables_bounds[i].lb;
      }
      else if (this->variables_bounds[i].ub < x[i]) {
         x[i] = this->variables_bounds[i].ub;
      }
   }
}

double Problem::compute_constraint_violation(double constraint, size_t j) const {
   return std::max(std::max(0., this->constraint_bounds[j].lb - constraint), constraint - this->constraint_bounds[j].ub);
}

// compute ||c||
double Problem::compute_constraint_violation(const std::vector<double>& constraints, Norm residual_norm) const {
   // create a lambda to avoid allocating an std::vector
   auto residual_function = [&](size_t j) {
      return this->compute_constraint_violation(constraints[j], j);
   };
   return norm(residual_function, constraints.size(), residual_norm);
}

// compute ||c_S|| for a given set of constraints
double Problem::compute_constraint_violation(const std::vector<double>& constraints, const std::vector<size_t>& constraint_set, Norm residual_norm)
const {
   auto residual_function = [&](size_t k) {
      size_t j = constraint_set[k];
      return this->compute_constraint_violation(constraints[j], j);
   };
   return norm(residual_function, constraint_set.size(), residual_norm);
}
bool Problem::is_constrained() const {
   return (0 < this->number_constraints);
}

/* native C++ problem */

//CppProblem::CppProblem(std::string name, int number_variables, int number_constraints, double (*objective)(std::vector<double> x), std::vector<double> (*objective_gradient)(std::vector<double> x)):
//Problem(name, number_variables, number_constraints),
//objective(objective),
//objective_gradient_(objective_gradient) {
//}
//
//double CppProblem::objective(std::vector<double>& x) {
//    return this->objective(x);
//}
//
//std::vector<double> CppProblem::objective_dense_gradient(std::vector<double>& x) {
//    return this->objective_gradient_(x);
//}
//
//SparseVector<double> CppProblem::objective_sparse_gradient(std::vector<double>& x) {
//    std::vector<double> dense_gradient = this->objective_gradient_(x);
//    SparseVector<double> sparse_gradient;
//    for (size_t i = 0; i < dense_gradient.size(); i++) {
//        if (dense_gradient[i] != 0.) {
//            sparse_gradient.insert(i, dense_gradient[i]);
//        }
//    }
//    return sparse_gradient;
//}
//
//double CppProblem::evaluate_constraint(int j, std::vector<double>& x) {
//    return this->constraints_[j](x);
//}
//
//std::vector<double> CppProblem::evaluate_constraints(std::vector<double>& x) {
//    std::vector<double> constraints(this->number_constraints);
//    for (int j = 0; j < this->number_constraints; j++) {
//        constraints[j] = this->evaluate_constraint(j, x);
//    }
//    return constraints;
//}
