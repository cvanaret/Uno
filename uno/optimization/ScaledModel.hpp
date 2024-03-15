// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALEDMODEL_H
#define UNO_SCALEDMODEL_H

#include "Model.hpp"
#include "preprocessing/Scaling.hpp"

class ScaledModel: public Model {
public:
   ScaledModel(std::unique_ptr<Model> constraint_index, Iterate& initial_iterate, const Options& options);

   [[nodiscard]] double get_variable_lower_bound(size_t variable_index) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t variable_index) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t constraint_index) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t constraint_index) const override;

   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] BoundType get_variable_bound_type(size_t variable_index) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t constraint_index) const override;
   [[nodiscard]] BoundType get_constraint_bound_type(size_t constraint_index) const override;

   [[nodiscard]] size_t get_number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t get_number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t get_number_hessian_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;
   void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override;

   [[nodiscard]] const std::vector<size_t>& get_linear_constraints() const override;

private:
   std::unique_ptr<Model> original_model;
   Scaling scaling;
};

inline ScaledModel::ScaledModel(std::unique_ptr<Model> original_model, Iterate& initial_iterate, const Options& options):
      Model(original_model->name + "_scaled", original_model->number_variables, original_model->number_constraints),
      original_model(std::move(original_model)),
      scaling(this->original_model->number_constraints, options.get_double("function_scaling_threshold")) {
   if (options.get_bool("scale_functions")) {
      // evaluate the gradients at the current point
      initial_iterate.evaluate_objective_gradient(*this->original_model);
      initial_iterate.evaluate_constraint_jacobian(*this->original_model);
      this->scaling.compute(initial_iterate.evaluations.objective_gradient, initial_iterate.evaluations.constraint_jacobian);
      // scale the gradients
      scale(initial_iterate.evaluations.objective_gradient, this->scaling.get_objective_scaling());
      for (size_t constraint_index: Range(this->original_model->number_constraints)) {
         scale(initial_iterate.evaluations.constraint_jacobian[constraint_index], this->scaling.get_constraint_scaling(constraint_index));
      }
   }
   // check the scaling factors
   assert(0 < this->scaling.get_objective_scaling() && "Objective scaling failed.");
   for ([[maybe_unused]] size_t constraint_index: Range(this->number_constraints)) {
      assert(0 < this->scaling.get_constraint_scaling(constraint_index) && "Constraint scaling failed.");
   }

   // the constraint repartition (inequality/equality, linear) is the same as in the original model
   this->equality_constraints.reserve(this->number_constraints);
   this->inequality_constraints.reserve(this->number_constraints);
   for (size_t constraint_index: this->original_model->equality_constraints) {
      this->equality_constraints.push_back(constraint_index);
   }
   for (size_t constraint_index: this->original_model->inequality_constraints) {
      this->inequality_constraints.push_back(constraint_index);
   }

   // the slacks are the same as in the original model
   this->original_model->slacks.for_each([&](size_t constraint_index, size_t variable_index) {
      this->slacks.insert(constraint_index, variable_index);
   });

   // the bounded variables are the same as in the original model
   for (size_t variable_index: this->original_model->lower_bounded_variables) {
      this->lower_bounded_variables.push_back(variable_index);
   }
   for (size_t variable_index: this->original_model->upper_bounded_variables) {
      this->upper_bounded_variables.push_back(variable_index);
   }
   for (size_t variable_index: this->original_model->single_lower_bounded_variables) {
      this->single_lower_bounded_variables.push_back(variable_index);
   }
   for (size_t variable_index: this->original_model->single_upper_bounded_variables) {
      this->single_upper_bounded_variables.push_back(variable_index);
   }
}

inline double ScaledModel::get_variable_lower_bound(size_t variable_index) const {
   return this->original_model->get_variable_lower_bound(variable_index);
}

inline double ScaledModel::get_variable_upper_bound(size_t variable_index) const {
   return this->original_model->get_variable_upper_bound(variable_index);
}

inline double ScaledModel::get_constraint_lower_bound(size_t constraint_index) const {
   const double lb = this->original_model->get_constraint_lower_bound(constraint_index);
   // scale
   return this->scaling.get_constraint_scaling(constraint_index)*lb;
}

inline double ScaledModel::get_constraint_upper_bound(size_t constraint_index) const {
   const double ub = this->original_model->get_constraint_upper_bound(constraint_index);
   // scale
   return this->scaling.get_constraint_scaling(constraint_index)*ub;
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
   for (size_t constraint_index: Range(this->number_constraints)) {
      constraints[constraint_index] *= this->scaling.get_constraint_scaling(constraint_index);
   }
}

inline void ScaledModel::evaluate_constraint_gradient(const std::vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const {
   this->original_model->evaluate_constraint_gradient(x, constraint_index, gradient);
   // scale
   scale(gradient, this->scaling.get_constraint_scaling(constraint_index));
}

inline void ScaledModel::evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
   this->original_model->evaluate_constraint_jacobian(x, constraint_jacobian);
   // scale
   for (size_t constraint_index: Range(this->number_constraints)) {
      scale(constraint_jacobian[constraint_index], this->scaling.get_constraint_scaling(constraint_index));
   }
}

inline void ScaledModel::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier,
      const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const {
   // scale the objective and constraint multipliers
   const double scaled_objective_multiplier = objective_multiplier*this->scaling.get_objective_scaling();
   // TODO preallocate this vector
   // TODO check if the multipliers should be scaled
   static std::vector<double> scaled_multipliers(this->number_constraints);
   for (size_t constraint_index: Range(this->number_constraints)) {
      scaled_multipliers[constraint_index] = scaling.get_constraint_scaling(constraint_index)*multipliers[constraint_index];
   }
   this->original_model->evaluate_lagrangian_hessian(x, scaled_objective_multiplier, scaled_multipliers, hessian);
}

inline BoundType ScaledModel::get_variable_bound_type(size_t variable_index) const {
   return this->original_model->get_variable_bound_type(variable_index);
}

inline FunctionType ScaledModel::get_constraint_type(size_t constraint_index) const {
   return this->original_model->get_constraint_type(constraint_index);
}

inline BoundType ScaledModel::get_constraint_bound_type(size_t constraint_index) const {
   return this->original_model->get_constraint_bound_type(constraint_index);
}

inline size_t ScaledModel::get_number_objective_gradient_nonzeros() const {
   return this->original_model->get_number_objective_gradient_nonzeros();
}

inline size_t ScaledModel::get_number_jacobian_nonzeros() const {
   return this->original_model->get_number_jacobian_nonzeros();
}

inline size_t ScaledModel::get_number_hessian_nonzeros() const {
   return this->original_model->get_number_hessian_nonzeros();
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
   for (size_t constraint_index: Range(iterate.number_constraints)) {
      iterate.multipliers.constraints[constraint_index] *= this->scaling.get_constraint_scaling(constraint_index) / this->scaling.get_objective_scaling();
   }

   // unscale the bound multipliers
   for (size_t variable_index: Range(iterate.number_variables)) {
      iterate.multipliers.lower_bounds[variable_index] /= this->scaling.get_objective_scaling();
      iterate.multipliers.upper_bounds[variable_index] /= this->scaling.get_objective_scaling();
   }
}

inline const std::vector<size_t>& ScaledModel::get_linear_constraints() const {
   return this->original_model->get_linear_constraints();
}

#endif // UNO_SCALEDMODEL_H
