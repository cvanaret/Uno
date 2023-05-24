// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BOUNDRELAXEDMODEL_H
#define UNO_BOUNDRELAXEDMODEL_H

#include "Model.hpp"

class BoundRelaxedModel: public Model {
public:
   BoundRelaxedModel(std::unique_ptr<Model> original_model, const Options& options);

   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;

   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] BoundType get_variable_bound_type(size_t i) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t j) const override;
   [[nodiscard]] BoundType get_constraint_bound_type(size_t j) const override;

   [[nodiscard]] size_t get_number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t get_number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t get_number_hessian_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;
   void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override;

   [[nodiscard]] const std::vector<size_t>& get_linear_constraints() const override;

private:
   std::unique_ptr<Model> original_model;
   const double relaxation_factor;
};

inline BoundRelaxedModel::BoundRelaxedModel(std::unique_ptr<Model> original_model, const Options& options):
      Model(original_model->name + "_boundrelaxed", original_model->number_variables, original_model->number_constraints),
      original_model(std::move(original_model)),
      relaxation_factor(options.get_double("tolerance")) {
   // the constraint repartition (inequality/equality, linear) is the same as in the original model
   this->equality_constraints.reserve(this->number_constraints);
   this->inequality_constraints.reserve(this->number_constraints);
   for (size_t j: this->original_model->equality_constraints) {
      this->equality_constraints.push_back(j);
   }
   for (size_t j: this->original_model->inequality_constraints) {
      this->inequality_constraints.push_back(j);
   }

   // the slacks are the same as in the original model
   this->original_model->slacks.for_each([&](size_t j, size_t i) {
      this->slacks.insert(j, i);
   });

   // the bounded variables are the same as in the original model
   for (size_t i: this->original_model->lower_bounded_variables) {
      this->lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->original_model->upper_bounded_variables) {
      this->upper_bounded_variables.push_back(i);
   }
   for (size_t i: this->original_model->single_lower_bounded_variables) {
      this->single_lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->original_model->single_upper_bounded_variables) {
      this->single_upper_bounded_variables.push_back(i);
   }
}

inline double BoundRelaxedModel::get_variable_lower_bound(size_t i) const {
   const double lower_bound = this->original_model->get_variable_lower_bound(i);
   // relax the bound
   return lower_bound - this->relaxation_factor * std::max(1., std::abs(lower_bound));
}

inline double BoundRelaxedModel::get_variable_upper_bound(size_t i) const {
   const double upper_bound = this->original_model->get_variable_upper_bound(i);
   // relax the bound
   return upper_bound + this->relaxation_factor * std::max(1., std::abs(upper_bound));
}

inline double BoundRelaxedModel::get_constraint_lower_bound(size_t j) const {
   return this->original_model->get_constraint_lower_bound(j);
}

inline double BoundRelaxedModel::get_constraint_upper_bound(size_t j) const {
   return this->original_model->get_constraint_upper_bound(j);
}

inline double BoundRelaxedModel::evaluate_objective(const std::vector<double>& x) const {
   return this->original_model->evaluate_objective(x);
}

inline void BoundRelaxedModel::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   this->original_model->evaluate_objective_gradient(x, gradient);
}

inline void BoundRelaxedModel::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_model->evaluate_constraints(x, constraints);
}

inline void BoundRelaxedModel::evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const {
   this->original_model->evaluate_constraint_gradient(x, j, gradient);
}

inline void BoundRelaxedModel::evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
   this->original_model->evaluate_constraint_jacobian(x, constraint_jacobian);
}

inline void BoundRelaxedModel::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier,
      const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const {
   this->original_model->evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
}

inline BoundType BoundRelaxedModel::get_variable_bound_type(size_t i) const {
   return this->original_model->get_variable_bound_type(i);
}

inline FunctionType BoundRelaxedModel::get_constraint_type(size_t j) const {
   return this->original_model->get_constraint_type(j);
}

inline BoundType BoundRelaxedModel::get_constraint_bound_type(size_t j) const {
   return this->original_model->get_constraint_bound_type(j);
}

inline size_t BoundRelaxedModel::get_number_objective_gradient_nonzeros() const {
   return this->original_model->get_number_objective_gradient_nonzeros();
}

inline size_t BoundRelaxedModel::get_number_jacobian_nonzeros() const {
   return this->original_model->get_number_jacobian_nonzeros();
}

inline size_t BoundRelaxedModel::get_number_hessian_nonzeros() const {
   return this->original_model->get_number_hessian_nonzeros();
}

inline void BoundRelaxedModel::get_initial_primal_point(std::vector<double>& x) const {
   this->original_model->get_initial_primal_point(x);
}

inline void BoundRelaxedModel::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_model->get_initial_dual_point(multipliers);
}

inline void BoundRelaxedModel::postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const {
   this->original_model->postprocess_solution(iterate, termination_status);
}

inline const std::vector<size_t>& BoundRelaxedModel::get_linear_constraints() const {
   return this->original_model->get_linear_constraints();
}

#endif // UNO_BOUNDRELAXEDMODEL_H
