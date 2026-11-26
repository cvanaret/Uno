// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "IQPSolver.hpp"
#include "LPSolver.hpp"
#include "QuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"

namespace uno {
   IQPSolver::IQPSolver(std::unique_ptr<LPSolver> qp_solver):
      SubproblemSolver(), qp_solver(std::move(qp_solver)) {
   }

   IQPSolver::~IQPSolver() = default;

   void IQPSolver::initialize_memory(const Subproblem& subproblem) {
      this->direction = Direction(subproblem.number_variables, subproblem.number_constraints);
      this->direction_SOC = Direction(subproblem.number_variables, subproblem.number_constraints);
      this->constraints_SOC.resize(subproblem.number_constraints);
      this->jacobian_vector_product_SOC.resize(subproblem.number_constraints);
      this->initial_point_SOC.resize(subproblem.number_variables);
      this->qp_solver->initialize_memory(subproblem);
   }

   void IQPSolver::compute_least_squares_multipliers(const Subproblem& /*subproblem*/, Iterate& /*iterate*/,
         Evaluations& /*evaluations*/, double /*multipliers_threshold*/) {
      DEBUG << "The IQP solver does not compute least-squares multipliers, keeping existing multipliers\n";
   }

   const Direction& IQPSolver::solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      this->direction.reset();
      // build the QuadraticProgram from the Subproblem at the current iterate
      QuadraticProgram& quadratic_program = this->qp_solver->get_quadratic_program();
      quadratic_program.fill(statistics, subproblem, current_iterate, trust_region_radius, current_evaluations, warmstart_information);

      // solve the QP
      this->qp_solver->solve(statistics, initial_point, this->direction, warmstart_information);
      // a solution with non-finite components (e.g. after a breakdown of the QP solver) cannot be used
      if (this->direction.status != SubproblemStatus::ERROR && !this->direction.is_finite()) {
         DEBUG << "The subproblem solution has non-finite components\n";
         this->direction.status = SubproblemStatus::ERROR;
      }

      // compute the dual direction
      compute_dual_direction(subproblem, current_iterate, this->direction.multipliers);
      direction.norm = norm_inf(view(direction.primals, 0, subproblem.problem.get_number_original_variables()));
      return this->direction;
   }

   bool IQPSolver::requires_hessian_matrix(const Subproblem& subproblem) const {
      return subproblem.has_curvature() && this->qp_solver->requires_hessian_matrix(subproblem);
   }

   // second-order corrections (Fletcher & Leyffer's filterSQP): the QP is solved again with the same objective, Jacobian
   // and trust region, the constants c(x_k) of the linearized constraints being replaced with c(x_k + d) - J_k d, where d
   // is the last (QP or SOC) step. The SOC step d_soc thus satisfies c(x_k + d) + J_k (d_soc - d) in the bounds
   bool IQPSolver::has_second_order_corrections() const {
      return true;
   }

   void IQPSolver::initialize_second_order_corrections(const Subproblem& subproblem, const Iterate& /*current_iterate*/,
         const Iterate& trial_iterate, Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      // the Jacobian J_k is that of the current iterate
      this->current_evaluations_SOC = &current_evaluations;
      this->compute_SOC_constraint_constants(subproblem, trial_iterate, this->direction, trial_evaluations);
   }

   const Direction& IQPSolver::compute_second_order_correction(Statistics& statistics, const Subproblem& subproblem,
         const Iterate& current_iterate) {
      this->direction_SOC.reset();
      // only the constraint bounds change (same active set estimate, Jacobian and trust region)
      QuadraticProgram& quadratic_program = this->qp_solver->get_quadratic_program();
      quadratic_program.set_constraint_constants(subproblem, this->constraints_SOC);
      WarmstartInformation warmstart_information{};
      warmstart_information.no_changes();
      warmstart_information.constraint_bounds_changed = true;
      this->initial_point_SOC.fill(0.);
      this->qp_solver->solve(statistics, this->initial_point_SOC, this->direction_SOC, warmstart_information);
      // restore the constants c(x_k): later solves at the same iterate (e.g. with a smaller radius) are unaffected
      quadratic_program.restore_constraint_constants(subproblem);
      if (this->direction_SOC.status != SubproblemStatus::ERROR && !this->direction_SOC.is_finite()) {
         DEBUG << "The SOC solution has non-finite components\n";
         this->direction_SOC.status = SubproblemStatus::ERROR;
      }
      compute_dual_direction(subproblem, current_iterate, this->direction_SOC.multipliers);
      this->direction_SOC.norm = norm_inf(view(this->direction_SOC.primals, 0, subproblem.problem.get_number_original_variables()));
      return this->direction_SOC;
   }

   void IQPSolver::update_second_order_corrections(const Subproblem& subproblem, const Iterate& trial_iterate,
         Evaluations& trial_evaluations) {
      // the next SOC linearizes around the last SOC step
      this->compute_SOC_constraint_constants(subproblem, trial_iterate, this->direction_SOC, trial_evaluations);
   }

   // c(x_k + d) - J_k d
   void IQPSolver::compute_SOC_constraint_constants(const Subproblem& subproblem, const Iterate& trial_iterate,
         const Direction& step, Evaluations& trial_evaluations) {
      assert(this->current_evaluations_SOC != nullptr && "The SOCs were not initialized");
      subproblem.problem.evaluate_constraints(trial_iterate, this->constraints_SOC.view(), trial_evaluations);
      this->jacobian_vector_product_SOC.fill(0.);
      subproblem.problem.add_jacobian_vector_product(view(step.primals, 0, subproblem.number_variables),
         this->jacobian_vector_product_SOC.view(), *this->current_evaluations_SOC);
      this->constraints_SOC -= this->jacobian_vector_product_SOC;
   }

   const SolverWorkspace& IQPSolver::get_workspace() const {
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
