// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIMALITYPROBLEM_H
#define UNO_OPTIMALITYPROBLEM_H

#include "OptimizationProblem.hpp"

class OptimalityProblem: public OptimizationProblem {
public:
   explicit OptimalityProblem(const Model& model);

   [[nodiscard]] double get_objective_multiplier() const override;
   void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const override;
   void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const override;
   void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const override;

   void set_infeasibility_measure(Iterate& iterate, Norm progress_norm) const override;
   void set_optimality_measure(Iterate& iterate) const override;
   [[nodiscard]] double compute_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Direction& direction,
         double step_length, Norm progress_norm) const override;
   [[nodiscard]] std::function<double(double)> compute_predicted_optimality_reduction_model(const Iterate& current_iterate,
         const Direction& direction, double step_length, const SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
   [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
   [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;

   [[nodiscard]] size_t number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t number_hessian_nonzeros() const override;
};

inline OptimalityProblem::OptimalityProblem(const Model& model):
      OptimizationProblem(model, model.number_variables, model.number_constraints) {
   // figure out bounded variables
   for (size_t variable_index: this->model.lower_bounded_variables) {
      this->lower_bounded_variables.push_back(variable_index);
   }
   for (size_t variable_index: this->model.upper_bounded_variables) {
      this->upper_bounded_variables.push_back(variable_index);
   }
   for (size_t variable_index: this->model.single_lower_bounded_variables) {
      this->single_lower_bounded_variables.push_back(variable_index);
   }
   for (size_t variable_index: this->model.single_upper_bounded_variables) {
      this->single_upper_bounded_variables.push_back(variable_index);
   }
}

inline double OptimalityProblem::get_objective_multiplier() const {
   return 1.;
}

inline void OptimalityProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
   iterate.evaluate_objective_gradient(this->model);
   objective_gradient = iterate.evaluations.objective_gradient;
}

inline void OptimalityProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
   iterate.evaluate_constraints(this->model);
   copy_from(constraints, iterate.evaluations.constraints);
}

inline void OptimalityProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
   iterate.evaluate_constraint_jacobian(this->model);
   constraint_jacobian = iterate.evaluations.constraint_jacobian;
}

inline void OptimalityProblem::evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   this->model.evaluate_lagrangian_hessian(x, this->get_objective_multiplier(), multipliers, hessian);
}

// infeasibility measure: constraint violation
inline void OptimalityProblem::set_infeasibility_measure(Iterate& iterate, Norm progress_norm) const {
   iterate.evaluate_constraints(this->model);
   iterate.progress.infeasibility = this->model.constraint_violation(iterate.evaluations.constraints, progress_norm);
}

// optimality measure: scaled objective
inline void OptimalityProblem::set_optimality_measure(Iterate& iterate) const {
   iterate.evaluate_objective(this->model);
   const double objective = iterate.evaluations.objective;
   iterate.progress.optimality = [=](double objective_multiplier) {
      return objective_multiplier*objective;
   };
}

inline double OptimalityProblem::compute_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Direction& direction,
      double step_length, Norm progress_norm) const {
   // predicted infeasibility reduction: "‖c(x)‖ - ‖c(x) + ∇c(x)^T (αd)‖"
   const double current_constraint_violation = this->model.constraint_violation(current_iterate.evaluations.constraints,
         progress_norm);
   const double trial_linearized_constraint_violation = this->model.linearized_constraint_violation(direction.primals,
         current_iterate.evaluations.constraints, current_iterate.evaluations.constraint_jacobian, step_length, progress_norm);
   return current_constraint_violation - trial_linearized_constraint_violation;
}

inline std::function<double(double)> OptimalityProblem::compute_predicted_optimality_reduction_model(const Iterate& current_iterate,
      const Direction& direction, double step_length, const SymmetricMatrix<double>& hessian) const {
   // predicted optimality reduction: "-∇f(x)^T (αd) - α^2/2 d^T H d"
   const double directional_derivative = dot(direction.primals, current_iterate.evaluations.objective_gradient);
   const double quadratic_product = hessian.quadratic_product(direction.primals, direction.primals);
   return [=](double objective_multiplier) {
      return step_length * (-objective_multiplier*directional_derivative) - step_length*step_length/2. * quadratic_product;
   };
}

inline double OptimalityProblem::variable_lower_bound(size_t variable_index) const {
   return this->model.variable_lower_bound(variable_index);
}

inline double OptimalityProblem::variable_upper_bound(size_t variable_index) const {
   return this->model.variable_upper_bound(variable_index);
}

inline double OptimalityProblem::constraint_lower_bound(size_t constraint_index) const {
   return this->model.constraint_lower_bound(constraint_index);
}

inline double OptimalityProblem::constraint_upper_bound(size_t constraint_index) const {
   return this->model.constraint_upper_bound(constraint_index);
}

inline size_t OptimalityProblem::number_objective_gradient_nonzeros() const {
   return this->model.get_number_objective_gradient_nonzeros();
}

inline size_t OptimalityProblem::number_jacobian_nonzeros() const {
   return this->model.get_number_jacobian_nonzeros();
}

inline size_t OptimalityProblem::number_hessian_nonzeros() const {
   return this->model.get_number_hessian_nonzeros();
}

#endif // UNO_OPTIMALITYPROBLEM_H
