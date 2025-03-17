// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LagrangeNewtonSubproblem.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
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

   void LagrangeNewtonSubproblem::evaluate_hessian(SymmetricMatrix<size_t, double>& hessian) {
      this->problem.evaluate_lagrangian_hessian(this->hessian_model, this->current_iterate.primals, this->current_multipliers.constraints, hessian);
   }

   void LagrangeNewtonSubproblem::evaluate_functions(SparseVector<double>& objective_gradient, Vector<double>& constraints,
         RectangularMatrix<double>& jacobian, SymmetricMatrix<size_t, double>& hessian, const WarmstartInformation& warmstart_information) {
      // objective gradient
      if (warmstart_information.objective_changed) {
         this->evaluate_objective_gradient(objective_gradient);
      }
      // constraints and Jacobian
      if (warmstart_information.constraints_changed) {
         this->evaluate_constraints(constraints);
      }
      if (warmstart_information.constraint_jacobian_changed) {
         this->evaluate_constraint_jacobian(jacobian);
      }
      // Lagrangian Hessian
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed || warmstart_information.constraint_jacobian_changed) {
         hessian.reset();
         this->evaluate_hessian(hessian);
      }
   }

   void LagrangeNewtonSubproblem::assemble_augmented_matrix(Statistics& statistics, SparseVector<double>& objective_gradient, Vector<double>& constraints,
         RectangularMatrix<double>& jacobian, SymmetricMatrix<size_t, double>& hessian, /* TODO remove */ SymmetricMatrix<size_t, double>& augmented_matrix,
         DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver, WarmstartInformation& warmstart_information) {
      this->evaluate_functions(objective_gradient, constraints, jacobian, hessian, warmstart_information);

      // TODO use matrix views
      if (warmstart_information.objective_changed || warmstart_information.constraint_jacobian_changed) {
         // form the KKT matrix
         augmented_matrix.set_dimension(this->number_variables + this->number_constraints);
         augmented_matrix.reset();
         // copy the Lagrangian Hessian in the top left block
         for (const auto [row_index, column_index, element]: hessian) {
            augmented_matrix.insert(element, row_index, column_index);
         }

         // Jacobian of general constraints
         for (size_t column_index: Range(this->number_constraints)) {
            for (const auto [row_index, derivative]: jacobian[column_index]) {
               augmented_matrix.insert(derivative, row_index, this->number_variables + column_index);
            }
            augmented_matrix.finalize_column(column_index);
         }
      }

      if (warmstart_information.hessian_sparsity_changed || warmstart_information.jacobian_sparsity_changed) {
         DEBUG << "Augmented matrix:\n" << augmented_matrix;
         DEBUG << "Performing symbolic analysis of the augmented matrix\n";
         linear_solver.do_symbolic_analysis(augmented_matrix);
         warmstart_information.hessian_sparsity_changed = warmstart_information.jacobian_sparsity_changed = false;
      }
      if (warmstart_information.objective_changed || warmstart_information.constraint_jacobian_changed) {
         DEBUG << "Performing numerical factorization of the augmented matrix\n";
         linear_solver.do_numerical_factorization(augmented_matrix);
         this->regularization_strategy.regularize_augmented_matrix(statistics, linear_solver, augmented_matrix, this->number_variables,
            this->number_constraints, this->problem.dual_regularization_parameter());
      }
   }

   void LagrangeNewtonSubproblem::assemble_augmented_rhs(SparseVector<double>& objective_gradient, Vector<double>& constraints,
         RectangularMatrix<double>& jacobian, Vector<double>& rhs, const WarmstartInformation& /*warmstart_information*/) const {
      // TODO use WarmstartInformation
      // Lagrangian gradient
      this->compute_lagrangian_gradient(objective_gradient, jacobian, rhs);
      for (size_t variable_index: Range(this->number_variables)) {
         rhs[variable_index] = -rhs[variable_index];
      }

      // constraints
      for (size_t constraint_index: Range(this->number_constraints)) {
         rhs[this->number_variables + constraint_index] = -constraints[constraint_index];
      }
      DEBUG2 << "RHS: "; print_vector(DEBUG2, view(rhs, 0, this->number_variables + this->number_constraints));
      DEBUG << '\n';
   }
} // namespace