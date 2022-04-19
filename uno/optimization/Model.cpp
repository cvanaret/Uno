#include <cmath>
#include <iostream>
#include <cassert>
#include <utility>
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"

bool is_finite(double value) {
   return std::abs(value) < std::numeric_limits<double>::infinity();
}

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
      linear_constraints(this->number_constraints) {
}

void Model::determine_bounds_types(std::vector<Range>& bounds, std::vector<ConstraintType>& status) {
   assert(bounds.size() == status.size());
   // build the "status" vector as a mapping (map/transform operation) of the "bounds" vector
   std::transform(begin(bounds), end(bounds), begin(status), [](const Range& bounds_i) {
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
      if (this->get_constraint_status(j) == EQUAL_BOUNDS) {
         this->equality_constraints.insert(j, current_equality_constraint);
         current_equality_constraint++;
      }
      else {
         this->inequality_constraints.insert(j, current_inequality_constraint);
         current_inequality_constraint++;
      }
   }
}

void Model::project_point_in_bounds(std::vector<double>& x) const {
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

// native C++ problem

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
