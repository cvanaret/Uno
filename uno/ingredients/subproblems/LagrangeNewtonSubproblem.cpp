// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LagrangeNewtonSubproblem.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   LagrangeNewtonSubproblem::LagrangeNewtonSubproblem(const OptimizationProblem& problem, Iterate& current_iterate,
      const Vector<double>& current_multipliers, HessianModel& hessian_model, double trust_region_radius):
            number_variables(problem.number_variables), number_constraints(problem.number_constraints),
            problem(problem), current_iterate(current_iterate), current_multipliers(current_multipliers), hessian_model(hessian_model),
            trust_region_radius(trust_region_radius) { }

   void LagrangeNewtonSubproblem::evaluate_objective_gradient(SparseVector<double>& objective_gradient) {
      this->problem.evaluate_objective_gradient(this->current_iterate, objective_gradient);
   }

   void LagrangeNewtonSubproblem::evaluate_constraints(Vector<double>& constraints) {
      this->problem.evaluate_constraints(this->current_iterate, constraints);
   }

   void LagrangeNewtonSubproblem::evaluate_constraint_jacobian(RectangularMatrix<double>& jacobian) {
      this->problem.evaluate_constraint_jacobian(this->current_iterate, jacobian);
   }

   void LagrangeNewtonSubproblem::evaluate_hessian(Statistics& statistics, SymmetricMatrix<size_t, double>& hessian) {
      this->hessian_model.evaluate(statistics, this->problem, this->current_iterate.primals, this->current_multipliers, hessian);
   }

   void LagrangeNewtonSubproblem::compute_lagrangian_gradient(SparseVector<double>& objective_gradient, RectangularMatrix<double>& jacobian,
         Vector<double>& gradient) const {
      gradient.fill(0.);

      // objective gradient
      for (const auto [variable_index, derivative]: objective_gradient) {
         gradient[variable_index] += derivative;
      }

      // constraints
      for (size_t constraint_index: Range(this->number_constraints)) {
         if (this->current_multipliers[constraint_index] != 0.) {
            for (const auto [variable_index, derivative]: jacobian[constraint_index]) {
               gradient[variable_index] -= this->current_multipliers[constraint_index] * derivative;
            }
         }
      }
   }

   double LagrangeNewtonSubproblem::dual_regularization_parameter() const {
      return this->problem.dual_regularization_parameter();
   }
} // namespace