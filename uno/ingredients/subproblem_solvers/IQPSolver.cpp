// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "IQPSolver.hpp"
#include "LPSolver.hpp"
#include "QuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   IQPSolver::IQPSolver(std::unique_ptr<LPSolver> qp_solver):
      SubproblemSolver(), qp_solver(std::move(qp_solver)) {
   }

   IQPSolver::~IQPSolver() = default;

   void IQPSolver::initialize_memory(const Subproblem& subproblem) {
      this->qp_solver->initialize_memory(subproblem);
   }

   void IQPSolver::solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& initial_point, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      // build the QuadraticProgram from the Subproblem at the current iterate
      QuadraticProgram& quadratic_program = this->qp_solver->get_quadratic_program();
      quadratic_program.fill(statistics, subproblem, trust_region_radius, current_evaluations, warmstart_information);

      // solve the QP
      this->qp_solver->solve(statistics, initial_point, direction, warmstart_information);

      // compute the dual direction
      compute_dual_direction(subproblem, direction.multipliers);
   }

   SolverWorkspace& IQPSolver::get_workspace() {
      return this->qp_solver->get_workspace();
   }

   // protected member functions

   // because of the way we form LPs/QPs, we get the new *multipliers* back from the solver. To get the dual
   // displacements/direction, we need to subtract the current multipliers
   void IQPSolver::compute_dual_direction(const Subproblem& subproblem, Multipliers& direction_multipliers) {
      view(direction_multipliers.constraints, 0, subproblem.number_constraints) -=
         view(subproblem.current_iterate.multipliers.constraints, 0, subproblem.number_constraints);
      view(direction_multipliers.lower_bounds, 0, subproblem.number_variables) -=
         view(subproblem.current_iterate.multipliers.lower_bounds, 0, subproblem.number_variables);
      view(direction_multipliers.upper_bounds, 0, subproblem.number_variables) -=
         view(subproblem.current_iterate.multipliers.upper_bounds, 0, subproblem.number_variables);
   }
} // namespace
