// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <stdexcept>
#include <utility>
#include "EQPSolver.hpp"
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "LinearSystem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "options/Options.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"

namespace uno {
   namespace {
      IterativeRefinementParameters read_iterative_refinement_parameters(const Options& options) {
         const IterativeRefinementParameters parameters{
            options.get_unsigned_int("iterative_refinement_min_steps"),
            options.get_unsigned_int("iterative_refinement_max_steps"),
            options.get_double("iterative_refinement_residual_ratio_max"),
            options.get_double("iterative_refinement_residual_ratio_singular"),
            options.get_double("iterative_refinement_improvement_factor")
         };
         if (parameters.max_steps < parameters.min_steps) {
            throw std::invalid_argument("The option iterative_refinement_max_steps must be larger than or equal to iterative_refinement_min_steps");
         }
         if (parameters.residual_ratio_singular < parameters.residual_ratio_max) {
            throw std::invalid_argument("The option iterative_refinement_residual_ratio_singular must not be smaller than iterative_refinement_residual_ratio_max");
         }
         return parameters;
      }
   } // namespace

   template <typename LowRankCorrection>
   EQPSolver<LowRankCorrection>::EQPSolver(const Options& options, LowRankCorrection correction):
         SubproblemSolver(),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options)),
         correction(std::move(correction)),
         iterative_refinement(read_iterative_refinement_parameters(options)) {
   }

   template <typename LowRankCorrection>
   EQPSolver<LowRankCorrection>::~EQPSolver() = default;

   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::initialize_memory(const Subproblem& subproblem) {
      if constexpr (LowRankCorrection::requires_hessian_matrix) {
         if (!subproblem.has_hessian_matrix()) {
            throw std::runtime_error("The subproblem does not have an explicit Hessian matrix and cannot be solved with a direct linear solver");
         }
      }
      this->direction = Direction(subproblem.number_variables, subproblem.number_constraints);
      // access the linear system of the linear solver
      auto& linear_system = this->linear_solver->get_linear_system();
      linear_system.initialize_augmented_system(subproblem);
      this->linear_solver->initialize_memory();
      // iterative refinement
      const size_t uncondensed_dimension = subproblem.problem.uncondensed_system_dimension();
      this->uncondensed_rhs.resize(uncondensed_dimension);
      this->uncondensed_solution.resize(uncondensed_dimension);
      this->uncondensed_residual.resize(uncondensed_dimension);
      this->uncondensed_correction.resize(uncondensed_dimension);
      this->augmented_rhs_buffer.resize(linear_system.dimension);
      this->augmented_solution_buffer.resize(linear_system.dimension);
   }

   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate, Evaluations& evaluations,
         double multipliers_threshold) {
      if (multipliers_threshold == 0.) {
         return;
      }
      DEBUG << "Computing least-squares multipliers\n";

      // compute least-square multipliers
      auto& linear_system = this->linear_solver->get_linear_system();

      // set up the linear system
      auto [hessian, primal_inertia_correction, jacobian, dual_inertia_correction] = subproblem.compute_block_augmented_matrix(linear_system);
      hessian.fill(0.); // no Hessian contribution
      primal_inertia_correction.fill(1.); // Identity block
      subproblem.evaluate_jacobian(iterate, jacobian, evaluations);
      dual_inertia_correction.fill(0.); // no dual regularization

      // perform the symbolic analysis once and for all
      if (!this->analysis_performed) {
         DEBUG << "Performing symbolic analysis of the indefinite system\n";
         this->linear_solver->do_symbolic_analysis();
         this->analysis_performed = true;
      }

      // factorize the matrix
      this->linear_solver->do_numerical_factorization(false);

      // assemble the RHS
      linear_system.rhs.fill(0.);
      evaluations.evaluate_objective_gradient(subproblem.problem.model, iterate.primals);
      view(linear_system.rhs.data(), subproblem.problem.model.number_variables) = evaluations.objective_gradient;
      for (size_t variable_index: Range(subproblem.problem.number_variables)) {
         linear_system.rhs[variable_index] -= (iterate.multipliers.lower_bounds[variable_index] +
            iterate.multipliers.upper_bounds[variable_index]);
      }

      // solve the linear system
      DEBUG3 << "KKT matrix values: " << linear_system.matrix_values << "\n";
      DEBUG3 << "RHS: " << linear_system.rhs << "\n\n";
      this->linear_solver->solve_indefinite_system(linear_system.solution.data());

      // set the constraint multipliers if their norm is reasonable
      const auto least_squares_multipliers = view(linear_system.solution.data(), subproblem.number_variables,
         subproblem.number_variables + subproblem.number_constraints);
      if (norm_inf(least_squares_multipliers) <= multipliers_threshold) {
         iterate.multipliers.constraints = least_squares_multipliers;
         DEBUG2 << "Least-squares multipliers set to " << least_squares_multipliers << '\n';
      }
      else {
         DEBUG << "Least-squares multipliers too large\n";
      }
   }

   template <typename LowRankCorrection>
   const Direction& EQPSolver<LowRankCorrection>::solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& /*initial_point*/, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("The direct linear solver does not support a trust region");
      }
      this->direction.reset();

      // access the linear system
      auto& linear_system = this->linear_solver->get_linear_system();

      // set up the linear system by evaluating the functions at the current iterate
      if (warmstart_information.new_iterate) {
         // assemble the augmented matrix
         auto block_augmented_matrix = subproblem.compute_block_augmented_matrix(linear_system);
         subproblem.evaluate_lagrangian_hessian(statistics, current_iterate, block_augmented_matrix.hessian);
         subproblem.evaluate_jacobian(current_iterate, block_augmented_matrix.jacobian, current_evaluations);

         // perform the symbolic analysis once and for all
         if (!this->analysis_performed) {
            DEBUG << "Performing symbolic analysis of the indefinite system\n";
            this->linear_solver->do_symbolic_analysis();
            this->analysis_performed = true;
         }

         // regularize the augmented matrix (this calls the analysis and the factorization)
         subproblem.regularize_augmented_matrix(statistics, subproblem.dual_regularization_factor(), *this->linear_solver,
            block_augmented_matrix);
         this->pivoting_tolerance_increased_for_current_matrix = false;
         if (!this->linear_solver->matrix_is_singular()) {
            this->correction.update(subproblem, *this->linear_solver);
         }

         // assemble the RHS
         subproblem.assemble_augmented_rhs(current_iterate, current_evaluations, linear_system.rhs);
         DEBUG3 << "RHS: " << linear_system.rhs << '\n';
      }

      // solve the linear system
      DEBUG3 << "KKT matrix values: " << linear_system.matrix_values << "\n\n";
      if (this->linear_solver->matrix_is_singular()) {
         this->linear_solver->solve_indefinite_system(linear_system.solution.data());
      }
      else {
         this->solve_with_iterative_refinement(&statistics, subproblem, current_iterate);
      }
      if (this->linear_solver->matrix_is_singular()) {
         this->direction.status = SubproblemStatus::INFEASIBLE;
      }
      else {
         this->correction.apply(linear_system.solution);
         // assemble the full primal-dual direction
         subproblem.assemble_primal_dual_direction(current_iterate, linear_system.solution, this->direction);
      }
      direction.norm = norm_inf(view(direction.primals, 0, subproblem.problem.get_number_original_variables()));
      return this->direction;
   }

   template <typename LowRankCorrection>
   bool EQPSolver<LowRankCorrection>::has_second_order_corrections() const {
      return true;
   }

   // equality constraints (lb == ub): residuals r(x) = c(x) - lb
   static void evaluate_constraint_residuals(const Subproblem& subproblem, const Iterate& iterate, Vector<double>& residuals,
         Evaluations& evaluations) {
      subproblem.problem.evaluate_constraints(iterate, residuals.view(), evaluations);
      residuals -= subproblem.problem.get_constraints_lower_bounds();
   }

   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::initialize_second_order_corrections(const Subproblem& subproblem, const Iterate& current_iterate,
         const Iterate& trial_iterate, Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      if (!this->SOC_initialized) {
         this->constraints_SOC.resize(subproblem.number_constraints);
         this->constraints_buffer_SOC.resize(subproblem.number_constraints);
         this->direction_SOC = Direction(subproblem.number_variables, subproblem.number_constraints);
         this->SOC_initialized = true;
      }
      // r_soc = r(x_trial) + α r(x_k)
      evaluate_constraint_residuals(subproblem, trial_iterate, this->constraints_SOC, trial_evaluations);
      evaluate_constraint_residuals(subproblem, current_iterate, this->constraints_buffer_SOC, current_evaluations);
      this->constraints_SOC += this->direction.primal_dual_step_length * this->constraints_buffer_SOC;
   }

   // precondition: the constraints have been evaluated at the trial iterate in trial_evaluations
   template <typename LowRankCorrection>
   const Direction& EQPSolver<LowRankCorrection>::compute_second_order_correction(const Subproblem& subproblem, const Iterate& current_iterate) {
      // access the linear system
      auto& linear_system = this->linear_solver->get_linear_system();

      // constraint part of the RHS: -r_soc (the primal part is reused from the last direction computation)
      auto rhs_constraints = view(linear_system.rhs, subproblem.number_variables, subproblem.number_variables +
         subproblem.number_constraints);
      rhs_constraints = -this->constraints_SOC;
      DEBUG3 << "SOC RHS: " << linear_system.rhs << '\n';

      // solve the linear system (the matrix is not modified: no singularity handling) and assemble the full primal-dual direction
      this->solve_with_iterative_refinement(nullptr, subproblem, current_iterate);
      this->correction.apply(linear_system.solution);
      subproblem.assemble_primal_dual_direction(current_iterate, linear_system.solution, this->direction_SOC);
      return this->direction_SOC;
   }

   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::update_second_order_corrections(const Subproblem& subproblem, const Iterate& trial_iterate,
         Evaluations& trial_evaluations) {
      // r_soc = α_soc r_soc + r(x_trial_soc)
      this->constraints_SOC.scale(this->direction_SOC.primal_dual_step_length);
      evaluate_constraint_residuals(subproblem, trial_iterate, this->constraints_buffer_SOC, trial_evaluations);
      this->constraints_SOC += this->constraints_buffer_SOC;
   }

   // solve the augmented system with the current factorization (RHS: linear_system.rhs, solution: linear_system.solution)
   // and refine the solution wrt the uncondensed Newton system, following IPOPT's PDFullSpaceSolver::Solve:
   // - refinement stops when the residual ratio is small enough (after at least min_steps steps), when max_steps steps
   //   were performed (max_steps = 0: no refinement step), or when the residual ratio stagnates;
   // - upon failure, the linear solver is asked to increase its pivoting tolerance (once per matrix) and the system is
   //   solved again;
   // - if that is not possible and the residual ratio is above residual_ratio_singular, the matrix is deemed singular and
   //   the inertia correction strategy increases the regularization (once per solve; requires statistics != nullptr).
   //   Otherwise, the current solution is accepted
   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::solve_with_iterative_refinement(Statistics* statistics, const Subproblem& subproblem,
         const Iterate& current_iterate) {
      auto& linear_system = this->linear_solver->get_linear_system();
      this->number_refinement_steps = 0;
      // the residuals require the full Hessian in the augmented matrix
      if constexpr (!LowRankCorrection::requires_hessian_matrix) {
         this->linear_solver->solve_indefinite_system(linear_system.solution.data());
         return;
      }
      const OptimizationProblem& problem = subproblem.problem;
      problem.uncondense_rhs(current_iterate, linear_system.rhs, this->uncondensed_rhs);

      bool matrix_assumed_singular = false;
      while (true) {
         // solve the augmented system and recover the uncondensed solution
         this->linear_solver->solve_indefinite_system(linear_system.rhs.data(), linear_system.solution.data(), 1);
         problem.expand_solution(current_iterate, this->uncondensed_rhs, linear_system.solution, this->uncondensed_solution);
         this->residual_ratio = this->compute_uncondensed_residual_ratio(subproblem, current_iterate);
         DEBUG2 << "Residual ratio: " << this->residual_ratio << '\n';

         // iterative refinement
         size_t number_steps = 0;
         bool refinement_succeeded = false;
         while (true) {
            if (this->iterative_refinement.min_steps <= number_steps && this->residual_ratio <= this->iterative_refinement.residual_ratio_max) {
               refinement_succeeded = true;
               break;
            }
            if (this->iterative_refinement.max_steps <= number_steps) {
               break;
            }
            this->compute_iterative_refinement_correction(subproblem, current_iterate);
            this->uncondensed_solution -= this->uncondensed_correction;
            ++number_steps;
            ++this->number_refinement_steps;
            const double previous_residual_ratio = this->residual_ratio;
            this->residual_ratio = this->compute_uncondensed_residual_ratio(subproblem, current_iterate);
            DEBUG2 << "Residual ratio after " << number_steps << " refinement step(s): " << this->residual_ratio << '\n';
            // stagnation
            if (this->iterative_refinement.min_steps < number_steps && this->iterative_refinement.residual_ratio_max < this->residual_ratio &&
                  this->iterative_refinement.improvement_factor * previous_residual_ratio < this->residual_ratio) {
               break;
            }
         }
         if (refinement_succeeded) {
            break;
         }
         DEBUG << "Iterative refinement failed with residual ratio " << this->residual_ratio << '\n';
         // the matrix was already deemed singular: live with the current solution
         if (matrix_assumed_singular) {
            break;
         }
         // first try to improve the quality of the factorization (once per matrix)
         if (!this->pivoting_tolerance_increased_for_current_matrix && this->linear_solver->increase_pivoting_tolerance()) {
            this->pivoting_tolerance_increased_for_current_matrix = true;
            DEBUG << "Refactorizing the augmented matrix with an increased pivoting tolerance\n";
            this->linear_solver->do_numerical_factorization(false);
            continue;
         }
         // otherwise, deem the matrix singular if the residuals are bad enough
         matrix_assumed_singular = true;
         if (this->residual_ratio < this->iterative_refinement.residual_ratio_singular) {
            DEBUG << "The current solution is accepted\n";
            break;
         }
         if (statistics == nullptr || !subproblem.regularize_singular_augmented_matrix(*statistics, subproblem.dual_regularization_factor(),
               *this->linear_solver, subproblem.compute_block_augmented_matrix(linear_system))) {
            DEBUG << "The matrix cannot be regularized further, the current solution is accepted\n";
            break;
         }
         // new matrix
         this->pivoting_tolerance_increased_for_current_matrix = false;
      }
      // write the refined solution of the augmented system
      view(linear_system.solution, 0, linear_system.dimension) = view(this->uncondensed_solution, 0, linear_system.dimension);
   }

   // ratio ||A u - r||_inf / (min(||u||_inf, 1e6 ||r||_inf) + ||r||_inf) of the uncondensed system A u = r (as in IPOPT)
   template <typename LowRankCorrection>
   double EQPSolver<LowRankCorrection>::compute_uncondensed_residual_ratio(const Subproblem& subproblem, const Iterate& current_iterate) {
      const auto& linear_system = this->linear_solver->get_linear_system();
      linear_system.compute_matrix_vector_product(this->uncondensed_solution.data(), this->uncondensed_residual.data());
      subproblem.problem.complete_uncondensed_matrix_vector_product(current_iterate, this->uncondensed_solution,
         this->uncondensed_residual);
      this->uncondensed_residual -= this->uncondensed_rhs;

      const double rhs_norm = norm_inf(this->uncondensed_rhs);
      const double solution_norm = norm_inf(this->uncondensed_solution);
      const double residual_norm = norm_inf(this->uncondensed_residual);
      DEBUG3 << "||rhs|| = " << rhs_norm << ", ||solution|| = " << solution_norm << ", ||residual|| = " << residual_norm << '\n';
      if (rhs_norm + solution_norm == 0.) {
         return residual_norm;
      }
      constexpr double max_condition_number = 1e6;
      return residual_norm / (std::min(solution_norm, max_condition_number * rhs_norm) + rhs_norm);
   }

   // solve A c = residual with the current factorization (the RHS of the linear system is preserved)
   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::compute_iterative_refinement_correction(const Subproblem& subproblem, const Iterate& current_iterate) {
      subproblem.problem.condense_rhs(current_iterate, this->uncondensed_residual, this->augmented_rhs_buffer);
      this->linear_solver->solve_indefinite_system(this->augmented_rhs_buffer.data(), this->augmented_solution_buffer.data(), 1);
      subproblem.problem.expand_solution(current_iterate, this->uncondensed_residual, this->augmented_solution_buffer,
         this->uncondensed_correction);
   }

   template <typename LowRankCorrection>
   const SolverWorkspace& EQPSolver<LowRankCorrection>::get_workspace() const {
      return this->linear_solver->get_linear_system();
   }

   template class EQPSolver<NoLowRankCorrection>;
   template class EQPSolver<WoodburyCorrection>;
} // namespace