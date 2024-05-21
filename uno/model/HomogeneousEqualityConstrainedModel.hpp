// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HOMOGENEOUSEQUALITYCONSTRAINEDMODEL_H
#define UNO_HOMOGENEOUSEQUALITYCONSTRAINEDMODEL_H

#include "Model.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/ChainCollection.hpp"
#include "symbolic/CollectionAdapter.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"

// generate an equality-constrained model by:
// - introducing slacks in inequality constraints
// - subtracting the (possibly nonzero) RHS of equality constraints
// all constraints are of the form "c(x) = 0"
class HomogeneousEqualityConstrainedModel: public Model {
public:
   explicit HomogeneousEqualityConstrainedModel(std::unique_ptr<Model> original_model);

   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override { return this->model->evaluate_objective(x); }
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override { this->model->evaluate_objective_gradient(x, gradient); }
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
   [[nodiscard]] BoundType get_variable_bound_type(size_t variable_index) const override;
   [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override { return this->lower_bounded_variables; }
   [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override { return this->upper_bounded_variables; }
   [[nodiscard]] const Collection<size_t>& get_slacks() const override { return this->slacks; }
   [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override { return this->single_lower_bounded_variables; }
   [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override { return this->single_upper_bounded_variables; }

   [[nodiscard]] double constraint_lower_bound(size_t /*constraint_index*/) const override { return 0.; } // c(x) = 0
   [[nodiscard]] double constraint_upper_bound(size_t /*constraint_index*/) const override { return 0.; }
   [[nodiscard]] FunctionType get_constraint_type(size_t constraint_index) const override { return this->model->get_constraint_type(constraint_index); }
   [[nodiscard]] BoundType get_constraint_bound_type(size_t /*constraint_index*/) const override { return EQUAL_BOUNDS; }
   [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override { return this->equality_constraints; }
   [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override { return this->inequality_constraints; }
   [[nodiscard]] const std::vector<size_t>& get_linear_constraints() const override { return this->model->get_linear_constraints(); }

   void initial_primal_point(std::vector<double>& x) const override;
   void initial_dual_point(std::vector<double>& multipliers) const override { this->model->initial_dual_point(multipliers); }
   void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override;

   [[nodiscard]] size_t number_objective_gradient_nonzeros() const override { return this->model->number_objective_gradient_nonzeros(); }
   [[nodiscard]] size_t number_jacobian_nonzeros() const override { return this->model->number_jacobian_nonzeros() + this->slacks.size(); }
   [[nodiscard]] size_t number_hessian_nonzeros() const override { return this->model->number_hessian_nonzeros(); }

protected:
   const std::unique_ptr<Model> model;
   std::vector<size_t> constraint_index_of_inequality_index;
   std::vector<size_t> slack_index_of_constraint_index;

   ForwardRange equality_constraints;
   ForwardRange inequality_constraints;
   SparseVector<size_t> slacks;
   std::vector<size_t> lower_bounded_slacks;
   std::vector<size_t> upper_bounded_slacks;
   std::vector<size_t> single_lower_bounded_slacks;
   std::vector<size_t> single_upper_bounded_slacks;
   ChainCollection<const Collection<size_t>&, CollectionAdapter<std::vector<size_t>&>> lower_bounded_variables;
   ChainCollection<const Collection<size_t>&, CollectionAdapter<std::vector<size_t>&>> upper_bounded_variables;
   ChainCollection<const Collection<size_t>&, CollectionAdapter<std::vector<size_t>&>> single_lower_bounded_variables;
   ChainCollection<const Collection<size_t>&, CollectionAdapter<std::vector<size_t>&>> single_upper_bounded_variables;
};

// Transform the problem into an equality-constrained problem with constraints c(x) = 0. This implies:
// - inequality constraints get a slack
// - equality constraints are shifted by their RHS
inline HomogeneousEqualityConstrainedModel::HomogeneousEqualityConstrainedModel(std::unique_ptr<Model> original_model):
      Model(original_model->name + "_equalityconstrained", original_model->number_variables + original_model->get_inequality_constraints().size(),
            original_model->number_constraints, original_model->objective_sign),
      // transfer ownership of the pointer
      model(std::move(original_model)),
      constraint_index_of_inequality_index(this->model->get_inequality_constraints().size()),
      slack_index_of_constraint_index(this->model->number_constraints),
      // all constraints are equality constraints
      equality_constraints(Range(this->number_constraints)),
      inequality_constraints(Range(0)),
      slacks(this->model->get_inequality_constraints().size()),
      lower_bounded_slacks(this->slacks.size()),
      upper_bounded_slacks(this->slacks.size()),
      lower_bounded_variables(concatenate(this->model->get_lower_bounded_variables(), CollectionAdapter(this->lower_bounded_slacks))),
      upper_bounded_variables(concatenate(this->model->get_upper_bounded_variables(), CollectionAdapter(this->upper_bounded_slacks))),
      single_lower_bounded_variables(concatenate(this->model->get_single_lower_bounded_variables(), CollectionAdapter(this->single_lower_bounded_slacks))),
      single_upper_bounded_variables(concatenate(this->model->get_single_upper_bounded_variables(), CollectionAdapter(this->single_upper_bounded_slacks))){
   // register the inequality constraint of each slack
   size_t inequality_index = 0;
   this->model->get_inequality_constraints().for_each([&](size_t, size_t constraint_index) {
      const size_t slack_variable_index = this->model->number_variables + inequality_index;
      this->constraint_index_of_inequality_index[inequality_index] = constraint_index;
      this->slack_index_of_constraint_index[constraint_index] = slack_variable_index;
      this->slacks.insert(constraint_index, slack_variable_index);
      if (is_finite(this->model->constraint_lower_bound(constraint_index))) {
         this->lower_bounded_slacks.push_back(slack_variable_index);
         if (not is_finite(this->model->constraint_upper_bound(constraint_index))) {
            this->single_lower_bounded_slacks.push_back(slack_variable_index);
         }
      }
      if (is_finite(this->model->constraint_upper_bound(constraint_index))) {
         this->upper_bounded_slacks.push_back(slack_variable_index);
         if (not is_finite(this->model->constraint_lower_bound(constraint_index))) {
            this->single_upper_bounded_slacks.push_back(slack_variable_index);
         }
      }
      inequality_index++;
   });
}

inline void HomogeneousEqualityConstrainedModel::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->model->evaluate_constraints(x, constraints);
   // inequality constraints: add the slacks
   this->get_slacks().for_each([&](size_t constraint_index, size_t slack_index) {
      constraints[constraint_index] -= x[slack_index];
   });

   // equality constraints: make sure they are homogeneous (c(x) = 0)
   this->model->get_equality_constraints().for_each([&](size_t, size_t constraint_index) {
      constraints[constraint_index] -= this->model->constraint_lower_bound(constraint_index);
   });
}

inline void HomogeneousEqualityConstrainedModel::evaluate_constraint_gradient(const std::vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const {
   this->model->evaluate_constraint_gradient(x, constraint_index, gradient);
   // if the original constraint is an inequality, add the slack contribution
   if (this->model->get_constraint_bound_type(constraint_index) != EQUAL_BOUNDS) {
      const size_t slack_variable_index = this->slack_index_of_constraint_index[constraint_index];
      gradient.insert(slack_variable_index, -1.);
   }
}

inline void HomogeneousEqualityConstrainedModel::evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
   this->model->evaluate_constraint_jacobian(x, constraint_jacobian);
   // add the slack contributions
   this->get_slacks().for_each([&](size_t constraint_index, size_t slack_index) {
      constraint_jacobian[constraint_index].insert(slack_index, -1.);
   });
}

inline void HomogeneousEqualityConstrainedModel::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   this->model->evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
   // extend the dimension of the Hessian by finalizing the remaining columns (note: the slacks do not enter the Hessian)
   for (size_t constraint_index: Range(this->model->number_variables, this->number_variables)) {
      hessian.finalize_column(constraint_index);
   }
}

inline double HomogeneousEqualityConstrainedModel::variable_lower_bound(size_t variable_index) const {
   if (variable_index < this->model->number_variables) { // original variable
      return this->model->variable_lower_bound(variable_index);
   }
   else { // slack variable
      const size_t slack_index = variable_index - this->model->number_variables;
      const size_t constraint_index = this->constraint_index_of_inequality_index[slack_index];
      return this->model->constraint_lower_bound(constraint_index);
   }
}

inline double HomogeneousEqualityConstrainedModel::variable_upper_bound(size_t variable_index) const {
   if (variable_index < this->model->number_variables) { // original variable
      return this->model->variable_upper_bound(variable_index);
   }
   else { // slack variable
      const size_t inequality_index = variable_index - this->model->number_variables;
      const size_t constraint_index = this->constraint_index_of_inequality_index[inequality_index];
      return this->model->constraint_upper_bound(constraint_index);
   }
}

inline BoundType HomogeneousEqualityConstrainedModel::get_variable_bound_type(size_t variable_index) const {
   if (variable_index < this->model->number_variables) { // original variable
      return this->model->get_variable_bound_type(variable_index);
   }
   else { // slack variable
      const size_t inequality_index = variable_index - this->model->number_variables;
      const size_t constraint_index = this->constraint_index_of_inequality_index[inequality_index];
      return this->model->get_constraint_bound_type(constraint_index);
   }
}

inline void HomogeneousEqualityConstrainedModel::initial_primal_point(std::vector<double>& x) const {
   this->model->initial_primal_point(x);
   // set the slacks
   this->get_slacks().for_each([&](size_t, size_t slack_index) {
      x[slack_index] = 0.;
   });
}

inline void HomogeneousEqualityConstrainedModel::postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const {
   this->model->postprocess_solution(iterate, termination_status);
   // discard the slacks
   iterate.number_variables = this->model->number_variables;
}

#endif // UNO_HOMOGENEOUSEQUALITYCONSTRAINEDMODEL_H
