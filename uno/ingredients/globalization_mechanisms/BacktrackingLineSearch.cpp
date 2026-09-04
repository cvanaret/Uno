// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "BacktrackingLineSearch.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/PredictedReductionModels.hpp"
#include "ingredients/inertia_correction_strategies/UnstableInertiaCorrection.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "optimization/EvaluationCache.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"
#include "tools/Symbols.hpp"

namespace uno {
   BacktrackingLineSearch::BacktrackingLineSearch(const Model& model, Options& options):
         GlobalizationMechanism(model, false, options),
         backtracking_ratio(options.get_double("LS_backtracking_ratio")),
         scale_duals_with_step_length(options.get_bool("LS_scale_duals_with_step_length")),
         delta(options.get_double("switching_delta")),
         gamma_alpha(options.get_double("gamma_alpha")),
         gamma_theta(1. - options.get_double("filter_beta")),
         gamma_phi(options.get_double("filter_gamma")),
         s_theta(options.get_double("switching_infeasibility_exponent")),
         s_phi(options.get_double("switching_objective_exponent")),
         theta_min(options.get_double("barrier_small_infeasibility_factor")),
         SOC_max_iterations(options.get_unsigned_int("SOC_max_iterations")),
         SOC_infeasibility_fraction(options.get_double("SOC_infeasibility_fraction")) {
      // check the initial and minimal step lengths
      assert(0 < this->backtracking_ratio && this->backtracking_ratio < 1. && "The LS backtracking ratio should be in (0, 1)");
   }

   void BacktrackingLineSearch::initialize(Statistics& statistics, const Model& model, Iterate& current_iterate,
         EvaluationCache& evaluation_cache, Options& options) {
      this->constraint_relaxation_strategy->initialize(statistics, current_iterate, false, evaluation_cache, options);
      statistics.add_column("LS", Statistics::int_width, 3);
      statistics.add_column("Steplength", Statistics::double_width + 1, 2);
      set_primal_statistics(statistics, model, current_iterate, evaluation_cache.current_evaluations);
      set_dual_residuals_statistics(statistics, current_iterate);
   }

   void BacktrackingLineSearch::compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, EvaluationCache& evaluation_cache, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) {
      DEBUG2 << "Current iterate\n" << current_iterate << '\n';

      // compute a feasible direction
      try {
         const Direction& direction = this->constraint_relaxation_strategy->compute_direction(statistics, current_iterate,
            INF<double>, evaluation_cache.current_evaluations, warmstart_information);
         if (direction.status != SubproblemStatus::INFEASIBLE) {
            check_unboundedness(direction);
            const PredictedReductionModels predicted_reduction_models =
               this->constraint_relaxation_strategy->build_predicted_reduction_models(current_iterate, direction,
               evaluation_cache.current_evaluations);
            const double minimum_step_length = compute_minimum_step_length(predicted_reduction_models, current_iterate, 1.);
            const bool backtracking_success = this->backtrack_along_direction(statistics, model, current_iterate, trial_iterate,
               direction, evaluation_cache, predicted_reduction_models, minimum_step_length, warmstart_information,
               user_callbacks);
            if (backtracking_success) {
               return;
            }
            // if backtracking failed, switch to the feasibility problem (below)
         }
         else { // infeasible
            statistics.set("Status", std::string("infeasible"));
            DEBUG << "/!\\ The subproblem is infeasible\n";
            // switch to the feasibility problem (below)
         }
      }
      // if the inertia correction failed, switch to solving the feasibility problem
      catch (const UnstableInertiaCorrection&) {
         statistics.set("Status", "inertia correction");
         // switch to the feasibility problem (below)
      }

      // feasibility problem: test corner cases
      if (model.number_constraints == 0) {
         throw std::runtime_error("The model is unconstrained, should not happen");
      }
      if (current_iterate.primal_infeasibility == 0.) {
         throw std::runtime_error("The current iterate is exactly feasible, should not happen");
      }
      // if the line search failed, switch to solving the feasibility problem (test first if we can)
      if (this->constraint_relaxation_strategy->solving_feasibility_problem()) {
         throw std::runtime_error("Line search failed");
      }

      // switch to solving the feasibility problem
      this->constraint_relaxation_strategy->switch_to_feasibility_problem(statistics, current_iterate,
         evaluation_cache.current_evaluations, warmstart_information);
      assert(this->constraint_relaxation_strategy->solving_feasibility_problem());
      const Direction& direction = this->constraint_relaxation_strategy->compute_direction(statistics, current_iterate,
         INF<double>, evaluation_cache.current_evaluations, warmstart_information);
      check_unboundedness(direction);
      const PredictedReductionModels predicted_reduction_models =
         this->constraint_relaxation_strategy->build_predicted_reduction_models(current_iterate, direction,
         evaluation_cache.current_evaluations);
      const double minimum_step_length = compute_minimum_step_length(predicted_reduction_models, current_iterate, 0.);
      const bool backtracking_success = this->backtrack_along_direction(statistics, model, current_iterate,
         trial_iterate, direction, evaluation_cache, predicted_reduction_models, minimum_step_length, warmstart_information,
         user_callbacks);
      if (!backtracking_success) {
         throw std::runtime_error("The line search failed");
      }
   }

   std::string BacktrackingLineSearch::get_name() const {
      return "LS " + this->constraint_relaxation_strategy->get_name();
   }

   // protected member functions

   void BacktrackingLineSearch::assemble_trial_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length) const {
      GlobalizationMechanism::assemble_trial_iterate(model, current_iterate, trial_iterate, direction,
         // primal step length
         step_length * direction.primal_dual_step_length,
         // constraint dual step length: scale or not with the LS step length
         (this->scale_duals_with_step_length ? step_length : 1.) * direction.primal_dual_step_length,
         // bound dual step length
         direction.bound_dual_step_length);
   }

   double BacktrackingLineSearch::compute_minimum_step_length(const PredictedReductionModels& predicted_reduction_models,
         const Iterate& current_iterate, double objective_multiplier) const {
      /*
      pred        = -∇φᵀd
      curr_theta = θ(x)            (constraint violation)

      alpha_min = γ_θ
      if pred > 0:
          alpha_min = min( γ_θ , γ_φ · θ / pred)
          if θ ≤ θ_min:
              alpha_min = min( alpha_min , δ · θ^{s_θ} / pred^{s_φ} )
      return γ_α · alpha_min
      */
      const double predicted_optimality_reduction = predicted_reduction_models.optimality_first_order_reduction(objective_multiplier);
      const double theta = current_iterate.progress.infeasibility;
      double alpha_min = this->gamma_theta;
      if (predicted_optimality_reduction > 0.) {
         alpha_min = std::min(alpha_min, this->gamma_phi * theta / predicted_optimality_reduction);
         if (theta <= this->theta_min) {
            alpha_min = std::min(alpha_min, this->delta * std::pow(theta, this->s_theta) /
               std::pow(predicted_optimality_reduction, this->s_phi));
         }
      }
      return this->gamma_alpha * alpha_min;
   }

   // go a fraction along the direction by finding an acceptable step length
   // returns true upon success, false upon failure
   bool BacktrackingLineSearch::backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, EvaluationCache& evaluation_cache,
         const PredictedReductionModels& predicted_reduction_models, double minimum_step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      double step_length = 1.;
      bool termination = false;
      size_t number_iterations = 0;
      DEBUG << "\nLine search: minimum step length set to " << minimum_step_length << '\n';
      while (!termination) {
         ++number_iterations;
         DEBUG << "\n\tLine-search iteration " << number_iterations << ", step_length " << step_length << '\n';
         if (1 < number_iterations) { statistics.start_new_line(); }
         statistics.set("Steplength", step_length);
         // the total step length is the step length scaled by the direction's step length
         const double total_step_length = step_length * direction.primal_dual_step_length;
         const ProgressMeasures predicted_reductions = predicted_reduction_models(total_step_length);

         bool is_acceptable = false;
         try {
            // take a step as a fraction of the direction
            assemble_trial_iterate(model, current_iterate, trial_iterate, direction, step_length);
            statistics.set("||Step||", step_length * direction.norm);

            is_acceptable = this->constraint_relaxation_strategy->is_iterate_acceptable(statistics, model, current_iterate,
               trial_iterate, direction, total_step_length, false, evaluation_cache.current_evaluations,
               evaluation_cache.trial_evaluations, predicted_reductions, warmstart_information, user_callbacks);
            set_primal_statistics(statistics, model, trial_iterate, evaluation_cache.trial_evaluations);

            // tiny direction test: if the primal direction is tiny over a certain number of successive iterations,
            // accept the step unconditionally
            if (!is_acceptable && number_iterations == 1 && is_tiny_direction(current_iterate, direction)) {
               // try to evaluate the functions at the trial iterate, so that the next subproblem is well defined
               evaluation_cache.trial_evaluations.evaluate_constraints(model, trial_iterate.primals);
               evaluation_cache.trial_evaluations.evaluate_jacobian(model, trial_iterate.primals);
               // TODO compute progress measures
               ++this->number_consecutive_tiny_directions;
               if (this->number_consecutive_tiny_directions >= this->consecutive_tiny_directions_threshold) {
                  is_acceptable = true;
                  DEBUG << "Accepting tiny step\n";
                  statistics.set("Status", std::string(symbols::check) + " (tiny)");
                  this->number_consecutive_tiny_directions = 0;
               }
            }
         }
         catch (const EvaluationError&) {
            statistics.set("Status", "eval. error");
            DEBUG << "An evaluation error occurred, reducing the step length further\n";
         }
         statistics.set("LS", number_iterations);

         // try second-order corrections if the full step was rejected
         if (!is_acceptable && number_iterations == 1 && this->constraint_relaxation_strategy->has_second_order_corrections() &&
               this->SOC_max_iterations >= 1 && trial_iterate.progress.infeasibility >= current_iterate.progress.infeasibility) {
            assert(step_length == 1.);
            is_acceptable = this->compute_second_order_directions(statistics, model, current_iterate, trial_iterate, direction,
               evaluation_cache, predicted_reductions, warmstart_information, user_callbacks);
         }

         if (is_acceptable) {
            termination = true;
         }
         // from here on, the trial iterate is rejected
         else {
            step_length = this->decrease_step_length(step_length);
            // note: if minimum_step_length = 0 and step_length reaches 0 by successive divisions, the strict inequality
            // protects from an endless loop
            if (step_length * direction.primal_dual_step_length > minimum_step_length) {
               // keep going
               evaluation_cache.trial_evaluations.reset();
            }
            else {
               // minimum step length reached
               DEBUG << "The line search failed with a step length <= " << minimum_step_length << '\n';
               // check if we can terminate at a first-order point
               if (trial_iterate.status != SolutionStatus::NOT_OPTIMAL) {
                  statistics.set("Status", "accepted (small step length)");
                  termination = true;
               }
               else {
                  // switch to solving the feasibility problem
                  statistics.set("Status", "small step length");
                  evaluation_cache.trial_evaluations.reset();
                  return false;
               }
            }
         }

         if (is_acceptable) {
            set_dual_residuals_statistics(statistics, trial_iterate);
         }
         if (Logger::level == INFO) statistics.print_current_line();
      } // end while loop
      return true;
   }

   bool BacktrackingLineSearch::is_tiny_direction(const Iterate& current_iterate, const Direction& direction) const {
      constexpr double macheps = std::numeric_limits<double>::epsilon();
      for (size_t variable_index: Range(current_iterate.number_variables)) {
         if (std::abs(direction.primals[variable_index]) / (1. + std::abs(current_iterate.primals[variable_index])) >= 10.*macheps) {
            return false;
         }
      }
      if (current_iterate.primal_infeasibility > this->theta_min) {
         return false;
      }
      return true;
   }

   bool BacktrackingLineSearch::compute_second_order_directions(Statistics& statistics, const Model& model,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, EvaluationCache& evaluation_cache,
         const ProgressMeasures& predicted_reductions, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) const {
      // enter second-order corrections
      DEBUG << "\nEntering second-order corrections\n";

      bool SOC_termination = false;
      bool is_acceptable = false;
      try {
         this->constraint_relaxation_strategy->initialize_second_order_corrections(current_iterate, trial_iterate,
            evaluation_cache.current_evaluations, evaluation_cache.trial_evaluations);
         double old_infeasibility_SOC = trial_iterate.progress.infeasibility;
         size_t SOC_iteration = 1;

         while (!SOC_termination) {
            DEBUG << "\n\tSOC iteration " << SOC_iteration << '\n';

            const Direction& direction_SOC = this->constraint_relaxation_strategy->compute_second_order_correction(current_iterate);
            assemble_trial_iterate(model, current_iterate, trial_iterate, direction_SOC, 1.);
            evaluation_cache.trial_evaluations.reset();

            is_acceptable = this->constraint_relaxation_strategy->is_iterate_acceptable(statistics, model, current_iterate,
               trial_iterate, direction /* this is correct, see IPOPT paper */, direction.primal_dual_step_length, false,
               evaluation_cache.current_evaluations, evaluation_cache.trial_evaluations, predicted_reductions,
               warmstart_information, user_callbacks);

            if (is_acceptable) {
               // terminate the SOCs and the backtracking
               SOC_termination = true;
               DEBUG << "SOC direction acceptable " << '\n';
            }
            else if (SOC_iteration >= this->SOC_max_iterations ||
                  trial_iterate.progress.infeasibility > this->SOC_infeasibility_fraction * old_infeasibility_SOC) {
               // terminate the SOCs and keep backtracking
               SOC_termination = true;
               DEBUG << "SOC done, resume backtracking " << '\n';
            }
            else {
               // continue the SOCs
               this->constraint_relaxation_strategy->update_second_order_corrections(trial_iterate, evaluation_cache.trial_evaluations);
               old_infeasibility_SOC = trial_iterate.progress.infeasibility;
               ++SOC_iteration;
               DEBUG << "SOC direction rejected, continue SOCs" << '\n';
            }
         }
      }
      catch (const EvaluationError&) {
         // terminate the SOCs and keep backtracking
         SOC_termination = true;
      }
      return is_acceptable;
   }

   // step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ...
   double BacktrackingLineSearch::decrease_step_length(double step_length) const {
      step_length *= this->backtracking_ratio;
      assert(0 < step_length && step_length <= 1 && "The line-search step length is not in (0, 1]");
      return step_length;
   }

   void BacktrackingLineSearch::check_unboundedness(const Direction& direction) {
      if (direction.status == SubproblemStatus::UNBOUNDED_PROBLEM) {
         throw std::runtime_error("The subproblem is unbounded, this should not happen. If the subproblem has curvature,"
            "use regularization. If not, use a trust-region method.\n");
      }
   }
} // namespace