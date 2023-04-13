// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALEDMODEL_H
#define UNO_SCALEDMODEL_H

#include "Model.hpp"
#include "preprocessing/Scaling.hpp"

class ScaledModel: public Model {
public:
   ScaledModel(std::unique_ptr<Model> original_model, Iterate& first_iterate, const Options& options);

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

   [[nodiscard]] size_t get_maximum_number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t get_maximum_number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t get_maximum_number_hessian_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;
   void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override;

private:
   std::unique_ptr<Model> original_model;
   Scaling scaling;
};

inline ScaledModel::ScaledModel(std::unique_ptr<Model> original_model, Iterate& first_iterate, const Options& options):
      Model(original_model->name + "_scaled", original_model->number_variables, original_model->number_constraints, original_model->problem_type),
      original_model(std::move(original_model)),
      scaling(this->original_model->number_constraints, options.get_double("function_scaling_threshold")) {
   if (options.get_bool("scale_functions")) {
      // evaluate the gradients at the current point
      first_iterate.evaluate_objective_gradient(*this->original_model);
      first_iterate.evaluate_constraint_jacobian(*this->original_model);
      this->scaling.compute(first_iterate.evaluations.objective_gradient, first_iterate.evaluations.constraint_jacobian);
      // scale the gradients
      scale(first_iterate.evaluations.objective_gradient, this->scaling.get_objective_scaling());
      for (size_t j: Range(this->original_model->number_constraints)) {
         scale(first_iterate.evaluations.constraint_jacobian[j], this->scaling.get_constraint_scaling(j));
      }
   }
   // check the scaling factors
   assert(0 < this->scaling.get_objective_scaling() && "Objective scaling failed.");
   for ([[maybe_unused]] size_t j: Range(this->number_constraints)) {
      assert(0 < this->scaling.get_constraint_scaling(j) && "Constraint scaling failed.");
   }

   // the constraint repartition (inequality/equality, linear) is the same as in the original model
   this->equality_constraints.reserve(this->number_constraints);
   this->inequality_constraints.reserve(this->number_constraints);
   this->linear_constraints.reserve(this->number_constraints);
   for (size_t j: this->original_model->equality_constraints) {
      this->equality_constraints.push_back(j);
   }
   for (size_t j: this->original_model->inequality_constraints) {
      this->inequality_constraints.push_back(j);
   }
   for (size_t j: this->original_model->linear_constraints) {
      this->linear_constraints.push_back(j);
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

inline double ScaledModel::get_variable_lower_bound(size_t i) const {
   return this->original_model->get_variable_lower_bound(i);
}

inline double ScaledModel::get_variable_upper_bound(size_t i) const {
   return this->original_model->get_variable_upper_bound(i);
}

inline double ScaledModel::get_constraint_lower_bound(size_t j) const {
   const double lb = this->original_model->get_constraint_lower_bound(j);
   // scale
   return this->scaling.get_constraint_scaling(j)*lb;
}

inline double ScaledModel::get_constraint_upper_bound(size_t j) const {
   const double ub = this->original_model->get_constraint_upper_bound(j);
   // scale
   return this->scaling.get_constraint_scaling(j)*ub;
}

inline double ScaledModel::evaluate_objective(const std::vector<double>& x) const {
   const double objective = this->original_model->evaluate_objective(x);
   // scale
   return this->scaling.get_objective_scaling()*objective;
}

inline void ScaledModel::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   this->original_model->evaluate_objective_gradient(x, gradient);
   // scale
   scale(gradient, this->scaling.get_objective_scaling());
}

inline void ScaledModel::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_model->evaluate_constraints(x, constraints);
   // scale
   for (size_t j: Range(this->number_constraints)) {
      constraints[j] *= this->scaling.get_constraint_scaling(j);
   }
}

inline void ScaledModel::evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const {
   this->original_model->evaluate_constraint_gradient(x, j, gradient);
   // scale
   scale(gradient, this->scaling.get_constraint_scaling(j));
}

inline void ScaledModel::evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
   this->original_model->evaluate_constraint_jacobian(x, constraint_jacobian);
   // scale
   for (size_t j: Range(this->number_constraints)) {
      scale(constraint_jacobian[j], this->scaling.get_constraint_scaling(j));
   }
}

inline void ScaledModel::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier,
      const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const {
   // scale the objective and constraint multipliers
   const double scaled_objective_multiplier = objective_multiplier*this->scaling.get_objective_scaling();
   // TODO preallocate this vector
   // TODO check if the multipliers should be scaled
   static std::vector<double> scaled_multipliers(this->number_constraints);
   for (size_t j: Range(this->number_constraints)) {
      scaled_multipliers[j] = scaling.get_constraint_scaling(j)*multipliers[j];
   }
   this->original_model->evaluate_lagrangian_hessian(x, scaled_objective_multiplier, scaled_multipliers, hessian);
}

inline BoundType ScaledModel::get_variable_bound_type(size_t i) const {
   return this->original_model->get_variable_bound_type(i);
}

inline FunctionType ScaledModel::get_constraint_type(size_t j) const {
   return this->original_model->get_constraint_type(j);
}

inline BoundType ScaledModel::get_constraint_bound_type(size_t j) const {
   return this->original_model->get_constraint_bound_type(j);
}

inline size_t ScaledModel::get_maximum_number_objective_gradient_nonzeros() const {
   return this->original_model->get_maximum_number_objective_gradient_nonzeros();
}

inline size_t ScaledModel::get_maximum_number_jacobian_nonzeros() const {
   return this->original_model->get_maximum_number_jacobian_nonzeros();
}

inline size_t ScaledModel::get_maximum_number_hessian_nonzeros() const {
   return this->original_model->get_maximum_number_hessian_nonzeros();
}

inline void ScaledModel::get_initial_primal_point(std::vector<double>& x) const {
   this->original_model->get_initial_primal_point(x);
}

inline void ScaledModel::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_model->get_initial_dual_point(multipliers);
}

inline void ScaledModel::postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const {
   this->original_model->postprocess_solution(iterate, termination_status);

   // unscale the objective value
   iterate.evaluations.objective /= this->scaling.get_objective_scaling();

   // unscale the constraint multipliers
   for (size_t j: Range(iterate.number_constraints)) {
      iterate.multipliers.constraints[j] *= this->scaling.get_constraint_scaling(j) / this->scaling.get_objective_scaling();
   }

   // unscale the bound multipliers
   for (size_t i: Range(iterate.number_variables)) {
      iterate.multipliers.lower_bounds[i] /= this->scaling.get_objective_scaling();
      iterate.multipliers.upper_bounds[i] /= this->scaling.get_objective_scaling();
   }
}

#endif // UNO_SCALEDMODEL_H
