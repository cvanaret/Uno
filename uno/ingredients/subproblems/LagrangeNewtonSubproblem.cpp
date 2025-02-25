// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LagrangeNewtonSubproblem.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   LagrangeNewtonSubproblem::LagrangeNewtonSubproblem(const OptimizationProblem& problem, Iterate& current_iterate,
      const Multipliers& current_multipliers, HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy,
      double trust_region_radius):
            number_variables(problem.number_variables), number_constraints(problem.number_constraints),
            problem(problem), current_iterate(current_iterate), current_multipliers(current_multipliers), hessian_model(hessian_model),
            regularization_strategy(regularization_strategy), trust_region_radius(trust_region_radius) { }

   void LagrangeNewtonSubproblem::evaluate_objective_gradient(SparseVector<double>& objective_gradient) {
      this->problem.evaluate_objective_gradient(this->current_iterate, objective_gradient);
   }

   void LagrangeNewtonSubproblem::evaluate_constraints(Vector<double>& constraints) {
      this->problem.evaluate_constraints(this->current_iterate, constraints);
   }

   void LagrangeNewtonSubproblem::evaluate_constraint_jacobian(RectangularMatrix<double>& jacobian) {
      this->problem.evaluate_constraint_jacobian(this->current_iterate, jacobian);
   }

   void LagrangeNewtonSubproblem::evaluate_hessian(SymmetricMatrix<size_t, double>& hessian) {
      this->hessian_model.evaluate(this->problem, this->current_iterate.primals, this->current_multipliers.constraints, hessian);
      // TODO evaluate with the following dependency:
      // Regularization(
      //	  InequalityHandlingMethod( # this possibly adds structured terms (diag barrier)
      //		HessianModel(
      //			ConstraintRelaxationStrategy(model)
      //		)
      //	)
      //)
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
         if (this->current_multipliers.constraints[constraint_index] != 0.) {
            for (const auto [variable_index, derivative]: jacobian[constraint_index]) {
               gradient[variable_index] -= this->current_multipliers.constraints[constraint_index] * derivative;
            }
         }
      }

      /*
      // bound multipliers
      for (size_t variable_index: Range(this->number_variables)) {
         gradient[variable_index] -= (this->current_multipliers.lower_bounds[variable_index] + this->current_multipliers.upper_bounds[variable_index]);
      }
       */
   }

   void LagrangeNewtonSubproblem::regularize_matrix(Statistics& statistics, DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver,
         SymmetricMatrix<size_t, double>& matrix) {
      this->regularization_strategy.regularize_matrix(statistics, linear_solver, matrix, this->number_variables, this->number_constraints,
         this->problem.dual_regularization_parameter());
   }
} // namespace