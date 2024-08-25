// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CUSTOM_MODEL_H
#define UNO_CUSTOM_MODEL_H

#include "Model.hpp"

template <typename Objective, typename ObjectiveGradient, typename Constraints, typename ConstraintJacobian, typename LagrangianHessian>
class CustomModel: public Model {
public:
   CustomModel(std::string name, size_t number_variables, size_t number_constraints, double objective_sign, const Objective& objective,
            const ObjectiveGradient& objective_gradient, const Constraints& constraints, const ConstraintJacobian& constraint_jacobian,
            const LagrangianHessian& lagrangian_hessian):
         Model(std::move(name), number_variables, number_constraints, objective_sign),
         custom_objective(objective), custom_objective_gradient(objective_gradient), custom_constraints(constraints),
         custom_constraint_jacobian(constraint_jacobian), custom_lagrangian_hessian(lagrangian_hessian) {
   }

   [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override {
      return this->custom_objective(x);
   }

   void evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const override {
      return this->custom_objective_gradient(x, gradient);
   }

   void evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const override {
      return this->custom_constraints(x, constraints);
   }

   void evaluate_constraint_gradient(const Vector<double>& /*x*/, size_t /*constraint_index*/, SparseVector<double>& /*gradient*/) const override {
      throw std::runtime_error("CustomModel::evaluate_constraint_gradient not implemented");
   }

   void evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override {
      return this->custom_constraint_jacobian(x, constraint_jacobian);
   }

   void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const override {
      return this->custom_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
   }

   // purely functions
   [[nodiscard]] double variable_lower_bound(size_t /*variable_index*/) const override {
      return -INF<double>;
   }

   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override {
      return variable_index == 0 ? 0.5 : INF<double>;
   }

   [[nodiscard]] BoundType get_variable_bound_type(size_t variable_index) const override {
      return variable_index == 0 ? BoundType::UNBOUNDED : BoundType::BOUNDED_UPPER;
   }

   [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override { return this->lower_bounded_variables; }
   [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override { return this->upper_bounded_variables; }
   [[nodiscard]] const SparseVector<size_t>& get_slacks() const override { return this->slacks; }
   [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override { return this->single_lower_bounded_variables; }
   [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override { return this->single_upper_bounded_variables; }

   [[nodiscard]] double constraint_lower_bound(size_t /*constraint_index*/) const override {
      return 0.;
   }

   [[nodiscard]] double constraint_upper_bound(size_t /*constraint_index*/) const override {
      return 0.;
   }

   [[nodiscard]] FunctionType get_constraint_type(size_t /*constraint_index*/) const override {
      return NONLINEAR;
   }

   [[nodiscard]] BoundType get_constraint_bound_type(size_t /*constraint_index*/) const override {
      return BoundType::BOUNDED_LOWER;
   }

   [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override { return this->equality_constraints; }
   [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override { return this->inequality_constraints; }
   [[nodiscard]] const std::vector<size_t>& get_linear_constraints() const override { return this->linear_constraints; }

   void initial_primal_point(Vector<double>& x) const override {
      x[0] = -2.;
      x[1] = 1.;
   }

   void initial_dual_point(Vector<double>& /*multipliers*/) const override {
   }

   void postprocess_solution(Iterate& /*iterate*/, TerminationStatus /*termination_status*/) const override {
   }

   [[nodiscard]] size_t number_objective_gradient_nonzeros() const override {
      return 0;
   }

   [[nodiscard]] size_t number_jacobian_nonzeros() const override {
      return 0;
   }
   [[nodiscard]] size_t number_hessian_nonzeros() const override {
      return 0;
   }

protected:
   const Objective* custom_objective;
   const ObjectiveGradient* custom_objective_gradient;
   const Constraints* custom_constraints;
   const ConstraintJacobian* custom_constraint_jacobian;
   const LagrangianHessian* custom_lagrangian_hessian;

   const ForwardRange lower_bounded_variables{0};
   const ForwardRange upper_bounded_variables{0};
   const SparseVector<size_t> slacks{};
   const ForwardRange single_lower_bounded_variables{0};
   const ForwardRange single_upper_bounded_variables{0};

   const ForwardRange equality_constraints{0};
   const ForwardRange inequality_constraints{0};
   const std::vector<size_t> linear_constraints{};
};

#endif // UNO_CUSTOM_MODEL_H