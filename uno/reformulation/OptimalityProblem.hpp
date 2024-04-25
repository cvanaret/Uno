// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIMALITYPROBLEM_H
#define UNO_OPTIMALITYPROBLEM_H

#include "OptimizationProblem.hpp"

class OptimalityProblem: public OptimizationProblem {
public:
   explicit OptimalityProblem(const Model& model);

   [[nodiscard]] double get_objective_multiplier() const override { return 1.; }
   void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const override;
   void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const override;
   void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] double variable_lower_bound(size_t variable_index) const override { return this->model.variable_lower_bound(variable_index); }
   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override { return this->model.variable_upper_bound(variable_index); }
   [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override { return this->model.get_lower_bounded_variables(); }
   [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override { return this->model.get_upper_bounded_variables(); }
   [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override { return this->model.get_single_lower_bounded_variables(); }
   [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override { return this->model.get_single_upper_bounded_variables(); }

   [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override { return this->model.constraint_lower_bound(constraint_index); }
   [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override { return this->model.constraint_upper_bound(constraint_index); }

   [[nodiscard]] size_t number_objective_gradient_nonzeros() const override { return this->model.number_objective_gradient_nonzeros(); }
   [[nodiscard]] size_t number_jacobian_nonzeros() const override { return this->model.number_jacobian_nonzeros(); }
   [[nodiscard]] size_t number_hessian_nonzeros() const override { return this->model.number_hessian_nonzeros(); }
};

inline OptimalityProblem::OptimalityProblem(const Model& model): OptimizationProblem(model, model.number_variables, model.number_constraints) {
}

inline void OptimalityProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
   iterate.evaluate_objective_gradient(this->model);
   // TODO change this
   objective_gradient = iterate.evaluations.objective_gradient;
}

inline void OptimalityProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
   iterate.evaluate_constraints(this->model);
   copy_from(constraints, iterate.evaluations.constraints);
}

inline void OptimalityProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
   iterate.evaluate_constraint_jacobian(this->model);
   // TODO change this
   constraint_jacobian = iterate.evaluations.constraint_jacobian;
}

inline void OptimalityProblem::evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   this->model.evaluate_lagrangian_hessian(x, this->get_objective_multiplier(), multipliers, hessian);
}

#endif // UNO_OPTIMALITYPROBLEM_H
