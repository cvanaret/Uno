// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "IQPSolver.hpp"
#include "LPSolver.hpp"
#include "QuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"

namespace uno {
   IQPSolver::IQPSolver(std::unique_ptr<LPSolver> qp_solver):
      SubproblemSolver(), qp_solver(std::move(qp_solver)) {
   }

   IQPSolver::~IQPSolver() = default;

   void IQPSolver::initialize_memory(const Subproblem& subproblem) {
      this->direction = Direction(subproblem.number_variables, subproblem.number_constraints);
      this->qp_solver->initialize_memory(subproblem);
   }

   void IQPSolver::generate_initial_iterate(const Subproblem& /*subproblem*/, Iterate& /*initial_iterate*/,
         Evaluations& /*evaluations*/) {
      // do nothing
   }

   void IQPSolver::compute_least_squares_multipliers(const Subproblem& /*subproblem*/, Iterate& /*iterate*/,
         Evaluations& /*evaluations*/) {
      DEBUG << "The EQP solver does not compute least-squares multipliers, keeping existing multipliers";
   }

   Direction& IQPSolver::solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      this->direction.reset();
      // build the QuadraticProgram from the Subproblem at the current iterate
      QuadraticProgram& quadratic_program = this->qp_solver->get_quadratic_program();
      quadratic_program.fill(statistics, subproblem, current_iterate, trust_region_radius, current_evaluations, warmstart_information);

      // solve the QP
      this->qp_solver->solve(statistics, initial_point, this->direction, warmstart_information);

      // compute the dual direction
      compute_dual_direction(subproblem, current_iterate, this->direction.multipliers);
      return this->direction;
   }

   bool IQPSolver::has_second_order_corrections() const {
      return false;
   }

   const Direction& IQPSolver::compute_second_order_correction(const Subproblem& /*subproblem*/, const Iterate& /*current_iterate*/,
         const Vector<double>& /*constraints_SOC*/) {
      throw std::runtime_error("No SOC implemented in IQPSolver");
   }

   SolverWorkspace& IQPSolver::get_workspace() const {
      return this->qp_solver->get_workspace();
   }

   // protected member functions

   // because of the way we form LPs/QPs, we get the new *multipliers* back from the solver. To get the dual
   // displacements/direction, we need to subtract the current multipliers
   void IQPSolver::compute_dual_direction(const Subproblem& subproblem, const Iterate& current_iterate, Multipliers& direction_multipliers) {
      view(direction_multipliers.constraints, 0, subproblem.number_constraints) -=
         view(current_iterate.multipliers.constraints, 0, subproblem.number_constraints);
      view(direction_multipliers.lower_bounds, 0, subproblem.number_variables) -=
         view(current_iterate.multipliers.lower_bounds, 0, subproblem.number_variables);
      view(direction_multipliers.upper_bounds, 0, subproblem.number_variables) -=
         view(current_iterate.multipliers.upper_bounds, 0, subproblem.number_variables);
   }
} // namespace
