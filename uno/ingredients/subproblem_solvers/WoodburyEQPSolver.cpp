// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "WoodburyEQPSolver.hpp"
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Direction.hpp"
#include "options/Options.hpp"

namespace uno {
   WoodburyEQPSolver::WoodburyEQPSolver(LBFGSHessian& hessian_model, const Options& options):
         SubproblemSolver(),
         hessian_model(hessian_model),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"))) {
   }

   void WoodburyEQPSolver::initialize_memory(const Subproblem& subproblem) {
      this->linear_solver->initialize_augmented_system(subproblem);
   }

   void WoodburyEQPSolver::solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& /*initial_point*/, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("The direct linear solver does not support a trust region");
      }

      // solve the subproblem with only the diagonal part
      // TODO build the subproblem from here
      this->linear_solver->solve_indefinite_system(statistics, subproblem, direction, current_evaluations, warmstart_information);
      if (this->linear_solver->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
      }

      // compute the low-rank corrections
   }

   SolverWorkspace& WoodburyEQPSolver::get_workspace() {
      return this->linear_solver->get_workspace();
   }
} // namespace