// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NONLINEARPROBLEM_H
#define UNO_NONLINEARPROBLEM_H

#include <string>
#include <vector>
#include <map>
#include "optimization/Iterate.hpp"
#include "optimization/Model.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "tools/Range.hpp"

class NonlinearProblem {
public:
   NonlinearProblem(const Model& model, size_t number_variables, size_t number_constraints);
   virtual ~NonlinearProblem() = default;

   const Model& model;
   const size_t number_variables; /*!< Number of variables */
   const size_t number_constraints; /*!< Number of constraints */

   [[nodiscard]] bool is_constrained() const;

   // SparseVector<size_t> equality_constraints{}; /*!< inequality constraints */
   std::vector<size_t> equality_constraints{}; /*!< inequality constraints */
   std::vector<size_t> inequality_constraints{}; /*!< inequality constraints */
   // lists of bounded variables
   std::vector<size_t> lower_bounded_variables{}; // indices of the lower-bounded variables
   std::vector<size_t> upper_bounded_variables{}; // indices of the upper-bounded variables
   std::vector<size_t> single_lower_bounded_variables{}; // indices of the single lower-bounded variables
   std::vector<size_t> single_upper_bounded_variables{}; // indices of the single upper-bounded variables

   // function evaluations
   [[nodiscard]] virtual double get_objective_multiplier() const = 0;
   [[nodiscard]] virtual double evaluate_objective(Iterate& iterate) const = 0;
   virtual void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const = 0;
   virtual void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const = 0;
   virtual void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const = 0;
   virtual void evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const = 0;

   [[nodiscard]] static double compute_optimality_stationarity_error(const Iterate& iterate, Norm residual_norm);
   [[nodiscard]] static double compute_feasibility_stationarity_error(const Iterate& iterate, Norm residual_norm);
   [[nodiscard]] double compute_complementarity_error(size_t number_variables, const std::vector<double>& primals,
         const std::vector<double>& constraints, const Multipliers& multipliers) const;
   [[nodiscard]] double compute_feasibility_complementarity_error(size_t number_variables, const std::vector<double>& primals,
         const std::vector<double>& constraints, const Multipliers& multipliers) const;

   [[nodiscard]] size_t get_number_original_variables() const;
   [[nodiscard]] virtual double get_variable_lower_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_variable_upper_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_constraint_lower_bound(size_t j) const = 0;
   [[nodiscard]] virtual double get_constraint_upper_bound(size_t j) const = 0;

   [[nodiscard]] virtual size_t get_maximum_number_objective_gradient_nonzeros() const = 0;
   [[nodiscard]] virtual size_t get_maximum_number_jacobian_nonzeros() const = 0;
   [[nodiscard]] virtual size_t get_maximum_number_hessian_nonzeros() const = 0;
};

inline NonlinearProblem::NonlinearProblem(const Model& model, size_t number_variables, size_t number_constraints):
      model(model), number_variables(number_variables), number_constraints(number_constraints) {
   this->equality_constraints.reserve(this->number_constraints);
   this->inequality_constraints.reserve(this->number_constraints);
}

inline bool NonlinearProblem::is_constrained() const {
   return (0 < this->number_constraints);
}

inline size_t NonlinearProblem::get_number_original_variables() const {
   return this->model.number_variables;
}

inline double NonlinearProblem::compute_optimality_stationarity_error(const Iterate& iterate, Norm residual_norm) {
   // norm of the Lagrangian gradient
   return norm(iterate.lagrangian_gradient, residual_norm);
}

inline double NonlinearProblem::compute_feasibility_stationarity_error(const Iterate& iterate, Norm residual_norm) {
   // norm of the constraints contribution of the Lagrangian gradient
   return norm(iterate.lagrangian_gradient.constraints_contribution, residual_norm);
}

// complementary slackness error
inline double NonlinearProblem::compute_complementarity_error(size_t number_variables, const std::vector<double>& primals,
      const std::vector<double>& constraints, const Multipliers& multipliers) const {
   double error = 0.;
   // bound constraints
   for (size_t i: Range(number_variables)) {
      if (0. < multipliers.lower_bounds[i]) {
         error = std::max(error, std::abs(multipliers.lower_bounds[i] * (primals[i] - this->get_variable_lower_bound(i))));
      }
      if (multipliers.upper_bounds[i] < 0.) {
         error = std::max(error, std::abs(multipliers.upper_bounds[i] * (primals[i] - this->get_variable_upper_bound(i))));
      }
   }
   // constraints
   for (size_t j: this->inequality_constraints) {
      if (0. < multipliers.constraints[j]) { // lower bound
         error = std::max(error, std::abs(multipliers.constraints[j] * (constraints[j] - this->get_constraint_lower_bound(j))));
      }
      else if (multipliers.constraints[j] < 0.) { // upper bound
         error = std::max(error, std::abs(multipliers.constraints[j] * (constraints[j] - this->get_constraint_upper_bound(j))));
      }
   }
   return error;
}

// complementary slackness error
inline double NonlinearProblem::compute_feasibility_complementarity_error(size_t number_variables, const std::vector<double>& primals,
      const std::vector<double>& constraints, const Multipliers& multipliers) const {
   double error = 0.;
   // bound constraints
   for (size_t i: Range(number_variables)) {
      if (0. < multipliers.lower_bounds[i]) {
         error = std::max(error, std::abs(multipliers.lower_bounds[i] * (primals[i] - this->get_variable_lower_bound(i))));
      }
      if (multipliers.upper_bounds[i] < 0.) {
         error = std::max(error, std::abs(multipliers.upper_bounds[i] * (primals[i] - this->get_variable_upper_bound(i))));
      }
   }
   // constraints
   for (size_t j: Range(constraints.size())) {
      // violated constraints
      if (constraints[j] < this->get_constraint_lower_bound(j)) { // lower violated
         error = std::max(error, std::abs((1. - multipliers.constraints[j]) * (constraints[j] - this->get_constraint_lower_bound(j))));
      }
      else if (this->get_constraint_upper_bound(j) < constraints[j]) { // upper violated
         error = std::max(error, std::abs((1. + multipliers.constraints[j]) * (constraints[j] - this->get_constraint_upper_bound(j))));
      }
      // satisfied constraints
      else if (0. < multipliers.constraints[j]) { // lower bound
         error = std::max(error, std::abs(multipliers.constraints[j] * (constraints[j] - this->get_constraint_lower_bound(j))));
      }
      else if (multipliers.constraints[j] < 0.) { // upper bound
         error = std::max(error, std::abs(multipliers.constraints[j] * (constraints[j] - this->get_constraint_upper_bound(j))));
      }
   }
   return error;
}

#endif // UNO_NONLINEARPROBLEM_H
