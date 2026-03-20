// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "EQPSolver.hpp"
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Direction.hpp"
#include "options/Options.hpp"

namespace uno {
   EQPSolver::EQPSolver(const Options& options):
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"))) {
   }

   void EQPSolver::initialize_memory(const Subproblem& subproblem) {
      this->linear_solver->initialize_augmented_system(subproblem);
   }

   void EQPSolver::solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& /*initial_point*/, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("The interior-point subproblem has a trust region. This is not implemented yet");
      }

      this->linear_solver->solve_indefinite_system(statistics, subproblem, direction, current_evaluations, warmstart_information);
      if (this->linear_solver->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
      }
   }

   SolverWorkspace& EQPSolver::get_workspace() {
      return this->linear_solver->get_workspace();
   }
} // namespace