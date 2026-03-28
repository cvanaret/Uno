// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "EQPSolver.hpp"
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "LinearSystem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   EQPSolver::EQPSolver(const Options& options):
         SubproblemSolver(),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"))) {
   }

   void EQPSolver::initialize_memory(const Subproblem& subproblem) {
      if (!subproblem.has_hessian_matrix()) {
         throw std::runtime_error("The subproblem does not have an explicit Hessian matrix and cannot be solved with a direct linear solver");
      }
      // access the linear system of the linear solver
      auto& linear_system = this->linear_solver->get_linear_system();
      linear_system.initialize_augmented_system(subproblem);
      this->linear_solver->initialize_memory();
   }

   void EQPSolver::solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& /*initial_point*/, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("The direct linear solver does not support a trust region");
      }

      // access the linear system
      auto& linear_system = this->linear_solver->get_linear_system();

      // set up the linear system by evaluating the functions at the current iterate
      if (warmstart_information.new_iterate) {
         // assemble the augmented matrix
         subproblem.evaluate_lagrangian_hessian(statistics, linear_system.matrix_values.data());
         const size_t number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
         subproblem.evaluate_jacobian(linear_system.matrix_values.data() + number_hessian_nonzeros, current_evaluations);

         // perform the symbolic analysis once and for all
         if (!this->analysis_performed) {
            DEBUG << "Performing symbolic analysis of the indefinite system\n";
            this->linear_solver->do_symbolic_analysis();
            this->analysis_performed = true;
         }

         // regularize the augmented matrix (this calls the analysis and the factorization)
         subproblem.regularize_augmented_matrix(statistics, linear_system.matrix_values.data(),
            subproblem.dual_regularization_factor(), *this->linear_solver);

         // assemble the RHS
         subproblem.assemble_augmented_rhs(current_evaluations, linear_system.rhs);
      }

      // solve the linear system
      this->linear_solver->solve_indefinite_system(linear_system.solution.data());
      if (this->linear_solver->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
         return;
      }
      // assemble the full primal-dual direction
      subproblem.assemble_primal_dual_direction(linear_system.solution, direction);
   }

   SolverWorkspace& EQPSolver::get_workspace() {
      return this->linear_solver->get_linear_system();
   }
} // namespace