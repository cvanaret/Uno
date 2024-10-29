// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FIXEDBOUNDSCONSTRAINTSMODEL_H
#define UNO_FIXEDBOUNDSCONSTRAINTSMODEL_H

#include "Model.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "symbolic/CollectionAdapter.hpp"
#include "symbolic/Concatenation.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   // move the fixed variables to the set of general constraints (required in barrier methods)
   // for instance, the constraint x[2] == 1 is not interpreted as 1 <= x[2] <= 1, but as the
   // linear constraint 0 x[1] + 1 x[2] + ... + 0 x[n] == 1
   class FixedBoundsConstraintsModel: public Model {
   public:
      FixedBoundsConstraintsModel(std::unique_ptr<Model> original_model, const Options& options);

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override { return this->model->evaluate_objective(x); }

      void evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const override {
         this->model->evaluate_objective_gradient(x, gradient);
      }

      void evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const override {
         this->model->evaluate_constraints(x, constraints);
         // add the fixed variables
         size_t current_constraint = this->model->number_constraints;
         for (size_t fixed_variable_index: this->model->get_fixed_variables()) {
            constraints[current_constraint] = x[fixed_variable_index];
            current_constraint++;
         }
      }

      void evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override {
         if (constraint_index < this->model->number_constraints) {
            // original constraint
            this->model->evaluate_constraint_gradient(x, constraint_index, gradient);
         }
         else {
            // fixed variable
            const size_t variable_index = this->model->get_fixed_variables()[constraint_index - this->model->number_constraints];
            gradient.insert(variable_index, 1.);
         }
      }

      void evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override {
         this->model->evaluate_constraint_jacobian(x, constraint_jacobian);
         // add the fixed variables
         size_t current_constraint = this->model->number_constraints;
         for (size_t fixed_variable_index: this->model->get_fixed_variables()) {
            constraint_jacobian[current_constraint].insert(fixed_variable_index, 1.);
            current_constraint++;
         }
      }

      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            SymmetricMatrix<size_t, double>& hessian) const override {
         this->model->evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
      }

      // only these two functions are redefined
      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override {
         if (this->model->variable_lower_bound(variable_index) == this->model->variable_upper_bound(variable_index)) {
            // remove bounds of fixed variables
            return -INF<double>;
         }
         return this->model->variable_lower_bound(variable_index);
      }

      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override {
         if (this->model->variable_lower_bound(variable_index) == this->model->variable_upper_bound(variable_index)) {
            // remove bounds of fixed variables
            return INF<double>;
         }
         return this->model->variable_upper_bound(variable_index);
      }

      [[nodiscard]] BoundType get_variable_bound_type(size_t variable_index) const override {
         if (this->model->variable_lower_bound(variable_index) == this->model->variable_upper_bound(variable_index)) {
            // fixed variable: remove the bounds
            return BoundType::UNBOUNDED;
         }
         else {
            return this->model->get_variable_bound_type(variable_index);
         }
      }

      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override {
         return this->lower_bounded_variables_collection;
      }

      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override {
         return this->upper_bounded_variables_collection;
      }

      [[nodiscard]] const SparseVector<size_t>& get_slacks() const override { return this->model->get_slacks(); }
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override { return this->model->get_single_lower_bounded_variables(); }
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override { return this->model->get_single_upper_bounded_variables(); }
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override { return this->fixed_variables; }

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override {
         if (constraint_index < this->model->number_constraints) {
            // original constraint
            return this->model->constraint_lower_bound(constraint_index);
         }
         else {
            // fixed variable
            const size_t variable_index = this->model->get_fixed_variables()[constraint_index - this->model->number_constraints];
            return this->model->variable_lower_bound(variable_index);
         }
      }

      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override {
         if (constraint_index < this->model->number_constraints) {
            // original constraint
            return this->model->constraint_upper_bound(constraint_index);
         }
         else {
            // fixed variable
            const size_t variable_index = this->model->get_fixed_variables()[constraint_index - this->model->number_constraints];
            return this->model->variable_lower_bound(variable_index);
         }
      }

      [[nodiscard]] FunctionType get_constraint_type(size_t constraint_index) const override {
         if (constraint_index < this->model->number_constraints) {
            // original constraint
            return this->model->get_constraint_type(constraint_index);
         }
         else {
            // fixed variables
            return FunctionType::LINEAR;
         }
      }

      [[nodiscard]] BoundType get_constraint_bound_type(size_t constraint_index) const override {
         if (constraint_index < this->model->number_constraints) {
            // original constraint
            return this->model->get_constraint_bound_type(constraint_index);
         }
         else {
            // fixed variables
            return BoundType::EQUAL_BOUNDS;
         }
      }

      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override { return this->equality_constraints; }
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override { return this->model->get_inequality_constraints(); }
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override { return this->linear_constraints; }

      void initial_primal_point(Vector<double>& x) const override {
         this->model->initial_primal_point(x);
         // set the fixed variables
         for (size_t variable_index: this->model->get_fixed_variables()) {
            x[variable_index] = this->model->variable_lower_bound(variable_index);
         }
      }
      void initial_dual_point(Vector<double>& multipliers) const override { this->model->initial_dual_point(multipliers); }
      void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override {
         this->model->postprocess_solution(iterate, termination_status);
      }

      [[nodiscard]] size_t number_objective_gradient_nonzeros() const override { return this->model->number_objective_gradient_nonzeros(); }
      [[nodiscard]] size_t number_jacobian_nonzeros() const override {
         return this->model->number_jacobian_nonzeros() + this->model->get_fixed_variables().size();
      }
      [[nodiscard]] size_t number_hessian_nonzeros() const override { return this->model->number_hessian_nonzeros(); }

   private:
      const std::unique_ptr<Model> model{};
      Vector<size_t> fixed_variables;
      Vector<size_t> lower_bounded_variables;
      CollectionAdapter<Vector<size_t>&> lower_bounded_variables_collection;
      Vector<size_t> upper_bounded_variables;
      CollectionAdapter<Vector<size_t>&> upper_bounded_variables_collection;
      Concatenation<const Collection<size_t>&, ForwardRange> equality_constraints;
      Concatenation<const Collection<size_t>&, ForwardRange> linear_constraints;
   };

   inline FixedBoundsConstraintsModel::FixedBoundsConstraintsModel(std::unique_ptr<Model> original_model, const Options& /*options*/):
         Model(original_model->name + "_fixedbounds", original_model->number_variables,
               // move the fixed variables to the set of general constraints
               original_model->number_constraints + original_model->get_fixed_variables().size(),
               original_model->objective_sign),
         model(std::move(original_model)),
         fixed_variables(), // empty vector
         lower_bounded_variables_collection(this->lower_bounded_variables),
         upper_bounded_variables_collection(this->upper_bounded_variables),
         equality_constraints(concatenate(this->model->get_equality_constraints(), Range(this->model->number_constraints, this->number_constraints))),
         linear_constraints(concatenate(this->model->get_linear_constraints(), Range(this->model->number_constraints, this->number_constraints))) {
      this->lower_bounded_variables.reserve(this->model->get_lower_bounded_variables().size());
      this->upper_bounded_variables.reserve(this->model->get_upper_bounded_variables().size());

      for (size_t variable_index: Range(this->model->number_variables)) {
         if (is_finite(this->variable_lower_bound(variable_index))) {
            this->lower_bounded_variables.emplace_back(variable_index);
         }
         if (is_finite(this->variable_upper_bound(variable_index))) {
            this->upper_bounded_variables.emplace_back(variable_index);
         }
      }
   }
} // namespace

#endif // UNO_FIXEDBOUNDSCONSTRAINTSMODEL_H