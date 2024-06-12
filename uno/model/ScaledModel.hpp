// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALEDMODEL_H
#define UNO_SCALEDMODEL_H

#include "Model.hpp"
#include "optimization/Iterate.hpp"
#include "preprocessing/Scaling.hpp"
#include "tools/Options.hpp"

class ScaledModel: public Model {
public:
   ScaledModel(std::unique_ptr<Model> original_model, Iterate& initial_iterate, const Options& options);

   [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
   void evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const override;
   void evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] double variable_lower_bound(size_t variable_index) const override { return this->model->variable_lower_bound(variable_index); }
   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override { return this->model->variable_upper_bound(variable_index); }
   [[nodiscard]] BoundType get_variable_bound_type(size_t variable_index) const override { return this->model->get_variable_bound_type(variable_index); }
   [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override { return this->model->get_lower_bounded_variables(); }
   [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override { return this->model->get_upper_bounded_variables(); }
   [[nodiscard]] const SparseVector<size_t>& get_slacks() const override { return this->model->get_slacks(); }
   [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override { return this->model->get_single_lower_bounded_variables(); }
   [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override { return this->model->get_single_upper_bounded_variables(); }

   [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
   [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t constraint_index) const override { return this->model->get_constraint_type(constraint_index); }
   [[nodiscard]] BoundType get_constraint_bound_type(size_t constraint_index) const override { return this->model->get_constraint_bound_type(constraint_index); }
   [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override { return this->model->get_equality_constraints(); }
   [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override { return this->model->get_inequality_constraints(); }
   [[nodiscard]] const std::vector<size_t>& get_linear_constraints() const override { return this->model->get_linear_constraints(); }

   void initial_primal_point(Vector<double>& x) const override { this->model->initial_primal_point(x); }
   void initial_dual_point(Vector<double>& multipliers) const override { this->model->initial_dual_point(multipliers); }
   void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override;

   [[nodiscard]] size_t number_objective_gradient_nonzeros() const override { return this->model->number_objective_gradient_nonzeros(); }
   [[nodiscard]] size_t number_jacobian_nonzeros() const override { return this->model->number_jacobian_nonzeros(); }
   [[nodiscard]] size_t number_hessian_nonzeros() const override { return this->model->number_hessian_nonzeros(); }

private:
   const std::unique_ptr<Model> model;
   Scaling scaling;
};

inline ScaledModel::ScaledModel(std::unique_ptr<Model> original_model, Iterate& initial_iterate, const Options& options):
      Model(original_model->name + "_scaled", original_model->number_variables, original_model->number_constraints, original_model->objective_sign),
      model(std::move(original_model)),
      scaling(this->model->number_constraints, options.get_double("function_scaling_threshold")) {
   if (options.get_bool("scale_functions")) {
      // evaluate the gradients at the current point
      initial_iterate.evaluate_objective_gradient(*this->model);
      initial_iterate.evaluate_constraint_jacobian(*this->model);
      this->scaling.compute(initial_iterate.evaluations.objective_gradient, initial_iterate.evaluations.constraint_jacobian);
      // scale the gradients
      scale(initial_iterate.evaluations.objective_gradient, this->scaling.get_objective_scaling());
      for (size_t constraint_index: Range(this->model->number_constraints)) {
         scale(initial_iterate.evaluations.constraint_jacobian[constraint_index], this->scaling.get_constraint_scaling(constraint_index));
      }
   }
   // check the scaling factors
   assert(0 < this->scaling.get_objective_scaling() && "Objective scaling failed.");
   for ([[maybe_unused]] size_t constraint_index: Range(this->number_constraints)) {
      assert(0 < this->scaling.get_constraint_scaling(constraint_index) && "Constraint scaling failed.");
   }
}

inline double ScaledModel::evaluate_objective(const Vector<double>& x) const {
   const double objective = this->model->evaluate_objective(x);
   return this->scaling.get_objective_scaling()*objective;
}

inline void ScaledModel::evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const {
   this->model->evaluate_objective_gradient(x, gradient);
   scale(gradient, this->scaling.get_objective_scaling());
}

inline void ScaledModel::evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const {
   this->model->evaluate_constraints(x, constraints);
   for (size_t constraint_index: Range(this->number_constraints)) {
      constraints[constraint_index] *= this->scaling.get_constraint_scaling(constraint_index);
   }
}

inline void ScaledModel::evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const {
   this->model->evaluate_constraint_gradient(x, constraint_index, gradient);
   scale(gradient, this->scaling.get_constraint_scaling(constraint_index));
}

inline void ScaledModel::evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
   this->model->evaluate_constraint_jacobian(x, constraint_jacobian);
   for (size_t constraint_index: Range(this->number_constraints)) {
      scale(constraint_jacobian[constraint_index], this->scaling.get_constraint_scaling(constraint_index));
   }
}

inline void ScaledModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   // scale the objective and constraint multipliers
   const double scaled_objective_multiplier = objective_multiplier*this->scaling.get_objective_scaling();
   // TODO preallocate this vector
   static Vector<double> scaled_multipliers(this->number_constraints);
   for (size_t constraint_index: Range(this->number_constraints)) {
      scaled_multipliers[constraint_index] = this->scaling.get_constraint_scaling(constraint_index) * multipliers[constraint_index];
   }
   this->model->evaluate_lagrangian_hessian(x, scaled_objective_multiplier, scaled_multipliers, hessian);
}

inline double ScaledModel::constraint_lower_bound(size_t constraint_index) const {
   const double lb = this->model->constraint_lower_bound(constraint_index);
   return this->scaling.get_constraint_scaling(constraint_index)*lb;
}

inline double ScaledModel::constraint_upper_bound(size_t constraint_index) const {
   const double ub = this->model->constraint_upper_bound(constraint_index);
   return this->scaling.get_constraint_scaling(constraint_index)*ub;
}

inline void ScaledModel::postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const {
   this->model->postprocess_solution(iterate, termination_status);

   // unscale the objective value
   if (iterate.is_objective_computed) {
      iterate.evaluations.objective /= this->scaling.get_objective_scaling();
   }

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

#endif // UNO_SCALEDMODEL_H
