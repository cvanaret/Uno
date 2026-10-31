// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include "EQPSolver.hpp"
#include "ingredients/subproblem/NewtonSystem.hpp"
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "LinearSystem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   namespace {
      IterativeRefinementParameters read_iterative_refinement_parameters(const Options& options) {
         const IterativeRefinementParameters parameters{
            options.get_unsigned_int("iterative_refinement_min_steps"),
            options.get_unsigned_int("iterative_refinement_max_steps"),
            options.get_double("iterative_refinement_residual_ratio_max"),
            options.get_double("iterative_refinement_residual_ratio_singular"),
            options.get_double("iterative_refinement_improvement_factor"),
         };
         // max_steps = 0 disables iterative refinement altogether (min_steps is then irrelevant)
         if (0 < parameters.max_steps && parameters.max_steps < parameters.min_steps) {
            throw std::invalid_argument("The option iterative_refinement_max_steps must be larger than or equal to iterative_refinement_min_steps");
         }
         if (!(0. <= parameters.residual_ratio_max)) {
            throw std::invalid_argument("The option iterative_refinement_residual_ratio_max must be nonnegative");
         }
         if (!(0. < parameters.improvement_factor && parameters.improvement_factor <= 1.)) {
            throw std::invalid_argument("The option iterative_refinement_improvement_factor must be in (0, 1]");
         }
         if (parameters.residual_ratio_singular < parameters.residual_ratio_max) {
            throw std::invalid_argument("The option iterative_refinement_residual_ratio_singular must not be smaller than iterative_refinement_residual_ratio_max");
         }
         return parameters;
      }
   } // namespace

   template <typename LowRankCorrection>
   EQPSolver<LowRankCorrection>::EQPSolver(const Options& options, LowRankCorrection low_rank_correction):
         SubproblemSolver(),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options)),
         low_rank_correction(std::move(low_rank_correction)),
         kkt_dumper(options),
         iterative_refinement(read_iterative_refinement_parameters(options)) {
      // if KKT dumper is enabled but the EQP solver does not use the explicit Hessian matrix, print a warning
      if constexpr (!LowRankCorrection::uses_explicit_hessian_matrix) {
         if (!options.get_string_optional("dump_kkt_path").value_or("").empty()) {
            WARNING << "'dump_kkt_path' is set, but KKT system dumping is only supported without low-rank corrections\n";
         }
      }
   }

   template <typename LowRankCorrection>
   EQPSolver<LowRankCorrection>::~EQPSolver() = default;

   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::initialize_memory(const Subproblem& subproblem) {
      if constexpr (LowRankCorrection::uses_explicit_hessian_matrix) {
         if (!subproblem.has_hessian_matrix()) {
            throw std::runtime_error("The subproblem does not have an explicit Hessian matrix and cannot be solved with a direct linear solver");
         }
      }
      this->direction = Direction(subproblem.number_variables, subproblem.number_constraints);
      this->zero_multipliers = Vector(subproblem.number_constraints, 0.);
      // symbolic Newton system and the linear system of the linear solver
      this->newton_system = subproblem.create_newton_system();
      auto& linear_system = this->linear_solver->get_linear_system();
      linear_system.initialize(*this->newton_system);
      this->linear_solver->initialize_memory();
      // iterative refinement on the uncondensed system
      const size_t uncondensed_dimension = this->newton_system->uncondensed_dimension();
      this->uncondensed_rhs.resize(uncondensed_dimension);
      this->uncondensed_solution.resize(uncondensed_dimension);
      this->uncondensed_residual.resize(uncondensed_dimension);
      this->uncondensed_correction.resize(uncondensed_dimension);
      this->condensed_rhs_buffer.resize(linear_system.dimension);
      this->condensed_solution_buffer.resize(linear_system.dimension);
   }

   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::initialize_statistics(Statistics& statistics) {
      if constexpr (LowRankCorrection::uses_explicit_hessian_matrix) {
         statistics.add_column("IR", Statistics::int_width, 0, /* is_extended = */ true);
      }
   }

   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::report_iterative_refinement(Statistics& statistics) {
      statistics.set("IR", std::to_string(this->iteration_refinement_steps) + join(this->iteration_refinement_flags, ","));
   }

   // least-squares multipliers: min ||∇L(x, y, z)|| wrt y, i.e. the solution of the system
   //    [I  J^T] [w]   [∇L(x, 0, z)]
   //    [J  0  ] [y] = [0          ]
   // where ∇L(x, 0, z) is the Lagrangian gradient of the subproblem's problem with zero constraint multipliers (∇f - z for
   // the barrier problem, whose barrier terms are represented by the bound multipliers z)
   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate, Evaluations& evaluations,
         double multipliers_threshold) {
      if (multipliers_threshold == 0.) {
         return;
      }
      DEBUG << "Computing least-squares multipliers\n";
      auto& linear_system = this->linear_solver->get_linear_system();

      // set up the linear system (this overwrites the augmented matrix and its factorization)
      this->factorization_holds_least_squares_system = true;
      auto [hessian, primal_inertia_correction, jacobian, dual_inertia_correction] =
         this->newton_system->block_augmented_matrix(linear_system.matrix_values.view());
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

      // factorize the matrix. If the Jacobian is rank deficient, regularize the (2, 2) block (when it is allocated)
      this->linear_solver->do_numerical_factorization(false);
      this->pivoting_tolerance_increased_for_current_matrix = false;
      if (this->linear_solver->matrix_is_singular() && !dual_inertia_correction.empty()) {
         DEBUG << "The least-squares system is singular, regularizing it\n";
         dual_inertia_correction.fill(-least_squares_dual_regularization);
         this->linear_solver->do_numerical_factorization(false);
      }
      // otherwise, the (consistent) singular system is solved anyway: a non-finite solution is rejected below

      // assemble the RHS: (∇L(x, 0, z), 0)
      auto rhs_lagrangian = view(linear_system.rhs, 0, subproblem.number_variables);
      subproblem.problem.evaluate_lagrangian_gradient(iterate.primals, this->zero_multipliers, iterate.multipliers.lower_bounds,
         iterate.multipliers.upper_bounds, evaluations, rhs_lagrangian);
      view(linear_system.rhs, subproblem.number_variables, linear_system.rhs.size()).fill(0.);

      // solve the linear system
      DEBUG3 << "KKT matrix values: " << linear_system.matrix_values << "\n";
      DEBUG3 << "RHS: " << linear_system.rhs << "\n\n";
      this->linear_solver->solve_indefinite_system(linear_system.solution.data());

      // set the constraint multipliers if their norm is reasonable (a NaN norm is rejected)
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
      this->iteration_refinement_steps = 0;
      this->iteration_refinement_flags.clear();

      // access the linear system
      auto& linear_system = this->linear_solver->get_linear_system();

      // set up the linear system by evaluating the functions at the current iterate
      if (warmstart_information.new_iterate || this->factorization_holds_least_squares_system) {
         this->factorization_holds_least_squares_system = false;
         // evaluate the Newton system (all the blocks except the regularization)
         auto block_augmented_matrix = this->newton_system->block_augmented_matrix(linear_system.matrix_values.view());
         this->newton_system->evaluate({statistics, current_iterate, current_evaluations}, linear_system.matrix_values.view());

         // perform the symbolic analysis once and for all
         if (!this->analysis_performed) {
            DEBUG << "Performing symbolic analysis of the indefinite system\n";
            this->linear_solver->do_symbolic_analysis();
            this->analysis_performed = true;
         }

         // regularize the augmented matrix (this calls the factorization)
         subproblem.regularize_augmented_matrix(statistics, subproblem.dual_regularization_factor(), *this->linear_solver,
            block_augmented_matrix);
         this->pivoting_tolerance_increased_for_current_matrix = false;
         if (!this->linear_solver->matrix_is_singular()) {
            this->low_rank_correction.update(subproblem, *this->linear_solver);
         }

         // assemble the RHS
         this->newton_system->assemble_rhs({statistics, current_iterate, current_evaluations}, linear_system.rhs.view());
         DEBUG3 << "RHS: " << linear_system.rhs << '\n';

         // optionally dump the regularized KKT system to disk
         if constexpr (LowRankCorrection::uses_explicit_hessian_matrix) {
            this->kkt_dumper.dump(linear_system, subproblem.problem.model.name,
               subproblem.problem.get_objective_multiplier() == 0., warmstart_information.iteration);
         }
      }

      // solve the linear system
      DEBUG3 << "KKT matrix values: " << linear_system.matrix_values << "\n\n";
      // solve with iterative refinement only if the matrix is not deemed singular
      if (this->linear_solver->matrix_is_singular()) {
         this->linear_solver->solve_indefinite_system(linear_system.solution.data());
      }
      else {
         this->solve_with_iterative_refinement(statistics, subproblem, /* allow_refactorization = */ true);
      }
      if (this->linear_solver->matrix_is_singular()) {
         this->direction.status = SubproblemStatus::INFEASIBLE;
      }
      else {
         this->low_rank_correction.apply(linear_system.solution);
         // assemble the full primal-dual direction
         subproblem.assemble_primal_dual_direction(current_iterate, *this->newton_system, linear_system.solution, this->direction);
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
   const Direction& EQPSolver<LowRankCorrection>::compute_second_order_correction(Statistics& statistics, const Subproblem& subproblem,
         const Iterate& current_iterate) {
      assert(!this->factorization_holds_least_squares_system && "The factorization is that of the least-squares system");
      // access the linear system
      auto& linear_system = this->linear_solver->get_linear_system();

      // constraint part of the RHS: -r_soc (the primal part is reused from the last direction computation)
      auto rhs_constraints = view(linear_system.rhs, subproblem.number_variables, subproblem.number_variables +
         subproblem.number_constraints);
      rhs_constraints = -this->constraints_SOC;
      DEBUG3 << "SOC RHS: " << linear_system.rhs << '\n';

      // solve the linear system (the factorization is not modified: no pivoting tolerance increase, no singularity
      // handling) and assemble the full primal-dual direction
      this->solve_with_iterative_refinement(statistics, subproblem, /* allow_refactorization = */ false);
      this->low_rank_correction.apply(linear_system.solution);
      subproblem.assemble_primal_dual_direction(current_iterate, *this->newton_system, linear_system.solution, this->direction_SOC);
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

   // solve the condensed system with the current factorization (RHS: linear_system.rhs, solution: linear_system.solution)
   // and refine the solution iteratively wrt the uncondensed Newton system, following IPOPT's PDFullSpaceSolver::Solve.
   // If the refined solution is not accurate:
   // 1. the pivoting tolerance of the linear solver is increased (at most once per matrix; it is never decreased) and the
   //    matrix is refactorized. The inertia correction strategy checks the inertia of the new factorization again and
   //    may increase the regularization;
   // 2. otherwise, if the residual ratio is below residual_ratio_singular, the solution is accepted;
   // 3. otherwise, the matrix is deemed singular (at most once per solve): the inertia correction strategy increases the
   //    regularization.
   // allow_refactorization == false: the factorization must not be modified (e.g. second-order corrections) and the
   // refined solution is accepted whatever its residual ratio.
   // max_steps = 0 (or a low-rank Hessian correction) disables the mechanism. Without refinement step, the solution is
   // that of the plain solve (bitwise)
   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::solve_with_iterative_refinement(Statistics& statistics, const Subproblem& subproblem,
         bool allow_refactorization) {
      auto& linear_system = this->linear_solver->get_linear_system();
      this->number_refinement_steps = 0;
      this->matrix_deemed_singular = false;
      // the residuals require the full Hessian in the condensed matrix (not available with low-rank corrections)
      if (!LowRankCorrection::uses_explicit_hessian_matrix || this->iterative_refinement.max_steps == 0) {
         this->residual_ratio = std::numeric_limits<double>::quiet_NaN(); // not computed
         this->linear_solver->solve_indefinite_system(linear_system.solution.data());
         return;
      }
      this->newton_system->uncondense_rhs(linear_system.rhs.view(), this->uncondensed_rhs.view());

      while (true) {
         // solve the condensed system, recover the uncondensed solution and refine it
         this->linear_solver->solve_indefinite_system(linear_system.solution.data());
         this->newton_system->expand_solution(this->uncondensed_rhs.view(), linear_system.solution.view(), this->uncondensed_solution.view());
         this->residual_ratio = this->compute_uncondensed_residual_ratio();
         DEBUG2 << "Residual ratio: " << this->residual_ratio << '\n';
         // 1. refinement reduced the residual below the threshold
         const bool refinement_success = this->refine_uncondensed_solution();
         if (refinement_success) {
            break;
         }
         DEBUG << "Iterative refinement failed with residual ratio " << this->residual_ratio << '\n';
         // the factorization cannot be modified, or the matrix was already deemed singular: live with the current solution
         if (!allow_refactorization || this->matrix_deemed_singular) {
            break;
         }
         const BlockAugmentedMatrix block_augmented_matrix = this->newton_system->block_augmented_matrix(linear_system.matrix_values.view());
         // 2. improve the quality of the factorization
         if (!this->pivoting_tolerance_increased_for_current_matrix && this->linear_solver->increase_pivoting_tolerance()) {
            this->pivoting_tolerance_increased_for_current_matrix = true;
            this->iteration_refinement_flags.emplace_back("piv");
            DEBUG << "Refactorizing the augmented matrix with an increased pivoting tolerance\n";
            subproblem.refactorize_augmented_matrix(statistics, subproblem.dual_regularization_factor(),
               *this->linear_solver, block_augmented_matrix);
            if (this->linear_solver->matrix_is_singular()) {
               // singular and not regularized: the caller handles the singularity
               DEBUG << "The refactorized augmented matrix is singular\n";
               this->linear_solver->solve_indefinite_system(linear_system.solution.data());
               this->report_iterative_refinement(statistics);
               return;
            }
            continue;
         }
         // 3. inaccurate, but acceptable solution
         if (this->residual_ratio < this->iterative_refinement.residual_ratio_singular) {
            DEBUG << "The current solution is accepted\n";
            this->iteration_refinement_flags.emplace_back("inac");
            break;
         }
         // 4. deem the matrix singular and regularize it
         this->matrix_deemed_singular = true;
         if (!subproblem.regularize_singular_augmented_matrix(statistics, subproblem.dual_regularization_factor(),
               *this->linear_solver, block_augmented_matrix)) {
            DEBUG << "The matrix cannot be regularized further, the current solution is accepted\n";
            break;
         }
         this->iteration_refinement_flags.emplace_back("sing");
      }
      this->report_iterative_refinement(statistics);
      // write the refined solution of the condensed system (the eliminated unknowns are recovered by back substitution)
      view(linear_system.solution, 0, linear_system.dimension) = view(this->uncondensed_solution, 0, linear_system.dimension);
   }

   // refinement steps on the current uncondensed solution (with an up-to-date residual ratio), as in IPOPT: stop when the
   // solution is accurate (after at least min_steps steps), after max_steps steps, or when the residual ratio stagnates
   // (after more than min_steps steps). Returns whether the solution is accurate
   template <typename LowRankCorrection>
   bool EQPSolver<LowRankCorrection>::refine_uncondensed_solution() {
      const IterativeRefinementParameters& parameters = this->iterative_refinement;
      size_t number_steps = 0;
      const auto accurate = [&] {
         return parameters.min_steps <= number_steps && this->residual_ratio <= parameters.residual_ratio_max;
      };
      while (!accurate() && number_steps < parameters.max_steps) {
         this->compute_iterative_refinement_correction();
         this->uncondensed_solution -= this->uncondensed_correction;
         ++number_steps;
         ++this->number_refinement_steps;
         ++this->iteration_refinement_steps;
         const double previous_residual_ratio = this->residual_ratio;
         this->residual_ratio = this->compute_uncondensed_residual_ratio();
         DEBUG2 << "Residual ratio after " << number_steps << " refinement step(s): " << this->residual_ratio << '\n';
         const bool stagnation = parameters.improvement_factor * previous_residual_ratio < this->residual_ratio;
         if (!accurate() && parameters.min_steps < number_steps && stagnation) {
            return false;
         }
      }
      return accurate();
   }

   // ratio ||K u - r||_inf / (min(||u||_inf, 1e6 ||r||_inf) + ||r||_inf) of the uncondensed system K u = r (as in IPOPT).
   // The product with K uses the factorized values (regularization included) of the kept block, without the Schur
   // complement (no cancellation), completed by the condensed blocks
   template <typename LowRankCorrection>
   double EQPSolver<LowRankCorrection>::compute_uncondensed_residual_ratio() {
      const auto& linear_system = this->linear_solver->get_linear_system();
      const auto [schur_begin, schur_end] = this->newton_system->schur_complement_nonzeros();
      linear_system.compute_matrix_vector_product(this->uncondensed_solution.data(), this->uncondensed_residual.data(), schur_begin,
         schur_end);
      this->newton_system->complete_uncondensed_product(this->uncondensed_solution.view(), this->uncondensed_residual.view());
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

   // solve K c = residual with the current factorization (the RHS of the linear system is preserved)
   template <typename LowRankCorrection>
   void EQPSolver<LowRankCorrection>::compute_iterative_refinement_correction() {
      this->newton_system->condense_rhs(this->uncondensed_residual.view(), this->condensed_rhs_buffer.view());
      this->linear_solver->solve_indefinite_system(this->condensed_rhs_buffer.data(), this->condensed_solution_buffer.data(), 1);
      this->newton_system->expand_solution(this->uncondensed_residual.view(), this->condensed_solution_buffer.view(),
         this->uncondensed_correction.view());
   }

   template <typename LowRankCorrection>
   const SolverWorkspace& EQPSolver<LowRankCorrection>::get_workspace() const {
      return this->linear_solver->get_linear_system();
   }

   template class EQPSolver<NoLowRankCorrection>;
   template class EQPSolver<WoodburyCorrection>;
} // namespace