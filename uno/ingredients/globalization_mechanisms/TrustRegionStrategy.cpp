// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <cassert>
#include <stdexcept>
#include "TrustRegionStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   TrustRegionStrategy::TrustRegionStrategy(const Model& model, const Options& options, std::vector<OptionOverride>& option_overrides) :
         GlobalizationMechanism(model, true, options, option_overrides),
         radius(options.get_double("TR_radius")),
         increase_factor(options.get_double("TR_increase_factor")),
         decrease_factor(options.get_double("TR_decrease_factor")),
         aggressive_decrease_factor(options.get_double("TR_aggressive_decrease_factor")),
         activity_tolerance(options.get_double("TR_activity_tolerance")),
         minimum_radius(options.get_double("TR_min_radius")),
         radius_reset_threshold(options.get_double("TR_radius_reset_threshold")),
         h_type_step_factor(options.get_double("TR_h_type_step_factor")),
         h_type_step_exponent(options.get_double("TR_h_type_step_exponent")),
         print_minor_iterations(options.get_bool("print_minor_iterations")) {
      if (!(0. < this->radius)) {
         throw std::runtime_error("The trust-region radius should be positive");
      }
      if (this->increase_factor <= 1.) {
         throw std::runtime_error("The trust-region increase factor should be > 1");
      }
      if (this->decrease_factor <= 1.) {
         throw std::runtime_error("The trust-region decrease factor should be > 1");
      }
      if (this->aggressive_decrease_factor <= 1.) {
         throw std::runtime_error("The trust-region aggressive decrease factor should be > 1");
      }
      if (!(0. < this->h_type_step_factor)) {
         throw std::runtime_error("The trust-region h-type step factor should be positive");
      }
      if (!(0.5 < this->h_type_step_exponent && this->h_type_step_exponent < 1.)) {
         throw std::runtime_error("The trust-region h-type step exponent should be in (0.5, 1)");
      }
      // the h-type step bound applies below the small-infeasibility threshold, which must be positive (h -> 0 along the
      // subsequence of FLT's Theorem 7 must eventually drop below it)
      if (!(0. < options.get_double("barrier_small_infeasibility_factor"))) {
         throw std::runtime_error("The small infeasibility factor should be positive with a trust region");
      }
   }

   void TrustRegionStrategy::initialize(Statistics& statistics, const Model& model, Iterate& current_iterate,
         EvaluationCache& evaluation_cache, const Options& options, std::vector<OptionOverride>& option_overrides) {
      this->constraint_relaxation_strategy->initialize(statistics, current_iterate, true, evaluation_cache, options,
         option_overrides);
      statistics.add_column("TR", Statistics::int_width, 3, /* is_extended = */ true);
      statistics.add_column("Radius", Statistics::double_width, 2, /* is_extended = */ true);
      statistics.set("Radius", this->radius);
      set_primal_statistics(statistics, model, current_iterate, evaluation_cache.current_evaluations);
      set_dual_residuals_statistics(statistics, current_iterate);
   }

   void TrustRegionStrategy::compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, EvaluationCache& evaluation_cache, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) {
      DEBUG2 << "Current iterate\n" << current_iterate << '\n';
      this->reset_radius();

      size_t number_iterations = 0;
      size_t number_consecutive_solver_errors = 0;
      bool termination = false;
      while (!termination) {
         bool is_acceptable = false;
         try {
            ++number_iterations;
            DEBUG << "\n\t### Trust-region inner iteration " << number_iterations << " with radius " << this->radius << "\n\n";
            if (this->print_minor_iterations && 1 < number_iterations && Logger::level == INFO) {
               statistics.start_new_line();
            }
            this->set_TR_statistics(statistics, number_iterations);

            // compute the direction within the trust region
            const Direction& direction = this->constraint_relaxation_strategy->compute_direction(statistics,
               current_iterate, this->radius, evaluation_cache.current_evaluations, warmstart_information);
            // a direction with non-finite components (e.g. after a breakdown of the subproblem solver) cannot be used
            SubproblemStatus direction_status = direction.status;
            if (direction_status != SubproblemStatus::INFEASIBLE && direction_status != SubproblemStatus::ERROR &&
                  !direction.is_finite()) {
               DEBUG << "The direction has non-finite components\n";
               direction_status = SubproblemStatus::ERROR;
            }
            if (direction_status == SubproblemStatus::INFEASIBLE) {
               if (model.number_constraints == 0) {
                  throw std::runtime_error("The model is unconstrained but the iterate is infeasible, should not happen");
               }
               if (this->constraint_relaxation_strategy->solving_feasibility_problem()) {
                  throw std::runtime_error("Uno is solving the feasibility problem but the subproblem is infeasible, should not happen");
               }
               statistics.set("Status", std::string("infeasible"));
               DEBUG << "/!\\ The subproblem is infeasible\n";

               // solve the feasibility problem at the next TR iteration with the same radius
               this->constraint_relaxation_strategy->switch_to_feasibility_problem(statistics, current_iterate,
                  evaluation_cache.current_evaluations, warmstart_information);
               continue;
            }
            // the subproblem solver failed repeatedly on the optimality subproblem at an infeasible iterate: the failure
            // does not depend on the radius. Solve the feasibility problem (a different subproblem) at the next TR
            // iteration with the same radius, as the line search does
            number_consecutive_solver_errors = (direction_status == SubproblemStatus::ERROR) ? number_consecutive_solver_errors + 1 : 0;
            if (TrustRegionStrategy::max_consecutive_solver_errors <= number_consecutive_solver_errors &&
                  0 < model.number_constraints && !this->constraint_relaxation_strategy->solving_feasibility_problem() &&
                  this->primal_tolerance < current_iterate.scaled_residuals.primal_infeasibility) {
               statistics.set("Status", "solver error");
               DEBUG << "/!\\ The subproblem solver failed repeatedly, switching to the feasibility problem\n";
               this->constraint_relaxation_strategy->switch_to_feasibility_problem(statistics, current_iterate,
                  evaluation_cache.current_evaluations, warmstart_information);
               number_consecutive_solver_errors = 0;
               continue;
            }
            assert(direction_status != SubproblemStatus::INFEASIBLE);
            statistics.set("||Step||", direction.norm);
            warmstart_information.no_changes();

            // deal with errors in the subproblem
            if (direction_status == SubproblemStatus::UNBOUNDED_PROBLEM) {
               // the subproblem is always bounded, but the objective may exceed a very large negative value
               statistics.set("Status", "unbounded subproblem");
               warmstart_information.trust_region_changed = true;
               this->decrease_radius_aggressively();
            }
            else if (direction_status == SubproblemStatus::ERROR) {
               statistics.set("Status", "solver error");
               this->decrease_radius();
               // reset the Hessian representation of the subproblem solver
               warmstart_information.whole_problem_changed();
            }
            else {
               // take full primal-dual step
               assemble_trial_iterate(current_iterate, trial_iterate, direction, direction.primal_dual_step_length,
                  direction.primal_dual_step_length, direction.bound_dual_step_length);
               evaluation_cache.trial_evaluations.reset();
               this->reset_active_trust_region_multipliers(model, direction, trial_iterate);

               const Decision decision = this->is_iterate_acceptable(statistics, model, current_iterate, trial_iterate,
                  direction, evaluation_cache, warmstart_information, user_callbacks);
               is_acceptable = decision.accepted;
               set_primal_statistics(statistics, model, trial_iterate, evaluation_cache.trial_evaluations);
               if (is_acceptable) {
                  // possibly increase the radius if trust region is active
                  this->possibly_increase_radius(direction.norm);
                  set_dual_residuals_statistics(statistics, trial_iterate);
                  termination = true;
               }
               else {
                  this->decrease_radius(direction.norm);
                  if (decision.reason == RejectionReason::H_TYPE_STEP_TOO_LONG) {
                     // geometric decrease with the bound b < ||d|| as a floor: the radii between b and the radius are
                     // sampled (jumping to b would skip FLT's window when the number of h-type iterations is finite), and
                     // the QP is solved again with a radius for which an h-type step is allowed
                     this->radius = std::max(this->radius, this->h_type_step_bound(current_iterate.progress.infeasibility));
                     DEBUG << "Trust-region radius set to " << this->radius << " (h-type step bound as a floor)\n";
                  }
                  warmstart_information.trust_region_changed = true;
                  // rejected step and small radius: try to terminate with a small step. The residuals of the rejected
                  // trial iterate were not computed: test (and terminate at) the current iterate instead, from which the
                  // trial iterate differs by a tiny step
                  if (this->radius < this->minimum_radius && this->check_termination_with_small_step(current_iterate)) {
                     trial_iterate = current_iterate;
                     evaluation_cache.trial_evaluations = evaluation_cache.current_evaluations;
                     is_acceptable = true;
                     termination = true;
                  }
               }
            }
         }
         // if an evaluation error occurs, decrease the radius
         catch (const EvaluationError&) {
            statistics.set("Status", "eval. error");
            DEBUG << "A function could not be evaluated. The trust-region radius will be reduced\n";
            this->decrease_radius();
            warmstart_information.trust_region_changed = true;
         }

         // print the statistics
         statistics.set("Time", Timer::format_to_seconds(statistics.timers.solve.get_elapsed_time()));
         if ((is_acceptable || this->print_minor_iterations) && Logger::level == INFO) {
            statistics.print_current_line();
         }

         // written so that a NaN radius also terminates the loop
         if (!is_acceptable && !(this->minimum_radius <= this->radius)) {
            throw std::runtime_error("Small radius");
         }
      }
   }

   void TrustRegionStrategy::postprocess_solution(Iterate& iterate, Evaluations& evaluations) const {
      this->constraint_relaxation_strategy->postprocess_solution(iterate, evaluations);
   }

   std::string TrustRegionStrategy::get_name() const {
      return "TR " + this->constraint_relaxation_strategy->get_name();
   }

   // protected member functions

   void TrustRegionStrategy::reset_active_trust_region_multipliers(const Model& model, const Direction& direction, Iterate& trial_iterate) const {
      assert(0 < this->radius && "The trust-region radius should be positive");
      // reset multipliers for bound constraints active at trust region (except if one of the original bounds is active)
      const auto& variables_lower_bounds = model.get_variables_lower_bounds();
      const auto& variables_upper_bounds = model.get_variables_upper_bounds();
      for (size_t variable_index: Range(model.number_variables)) {
         if (std::abs(direction.primals[variable_index] + this->radius) <= this->activity_tolerance &&
               this->activity_tolerance < std::abs(trial_iterate.primals[variable_index] - variables_lower_bounds[variable_index])) {
            trial_iterate.multipliers.lower_bounds[variable_index] = 0.;
         }
         if (std::abs(direction.primals[variable_index] - this->radius) <= this->activity_tolerance &&
               this->activity_tolerance < std::abs(variables_upper_bounds[variable_index] - trial_iterate.primals[variable_index])) {
            trial_iterate.multipliers.upper_bounds[variable_index] = 0.;
         }
      }
   }

   // h-type step bound b = factor h^nu, applied to the optimality phase if h <= theta (small-infeasibility threshold).
   // FLT (Theorem 7, h-type case) relies on the monotonicity of the predicted reduction in the radius, a consequence of the
   // global optimality of the QP solution. With KKT points of nonconvex QPs, a long h-type step at a large radius could be
   // accepted near a feasible non-KKT point, infinitely often. Any KKT point of the QP is f-type and acceptable for radii
   // in FLT's window [mu h, rho_bar] with rho_bar = min(sqrt(2 beta tau / (mnM)), kappa) (Lemma 5, with
   // q(d) <= q(d_hat) + 2nM rho^2 instead of global optimality). Since mu h << h^nu << sqrt(h) for nu in (1/2, 1) and small
   // h, the radii [b, 2b] with b = factor h^nu lie in the window: forbidding h-type steps longer than b lets the inner loop
   // sample an f-type step first. A KKT point d with ||d||_inf <= b is also a KKT point for the radius b (the trust-region
   // constraints are inactive), hence the test on the step norm rather than the radius.
   // After such a rejection, the radius decreases geometrically with b as a floor. A rejected step at a radius larger than
   // rho_bar has ||d|| > rho_bar (otherwise, d would be an acceptable f-type KKT point for a radius in the window), so the
   // next radius exceeds rho_bar/2: the accepted radii stay bounded away from 0, as required when the number of h-type
   // iterations is finite. Jumping straight to b (-> 0) would skip the window.
   // The threshold theta > 0 is fixed: along the subsequence of Theorem 7, h -> 0 eventually drops below it, and the case
   // of finitely many h-type iterations only relies on the radius update. Likewise, the proof only uses b for h -> 0, so
   // the factor is arbitrary: theta and factor = 100 avoid shortening h-type steps far from feasibility (at h = theta,
   // b = 0.1), where they are legitimate
   double TrustRegionStrategy::h_type_step_bound(double infeasibility) const {
      return std::max(this->h_type_step_factor * std::pow(infeasibility, this->h_type_step_exponent), this->minimum_radius);
   }

   // the trial iterate is accepted by the constraint relaxation strategy or if the step is small and we cannot switch to solving the feasibility problem
   Decision TrustRegionStrategy::is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, EvaluationCache& evaluation_cache,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      const PredictedReductionModels predicted_reduction_models =
         this->constraint_relaxation_strategy->build_predicted_reduction_models(current_iterate, direction,
         evaluation_cache.current_evaluations);
      const ProgressMeasures predicted_reductions = predicted_reduction_models(/* step_length = */ 1.,
         /* use_curvature_information = */ true);
      // the bound only applies to the optimality phase (FLT treat restoration as a black box), below the small-infeasibility
      // threshold theta > 0 of that phase
      const double infeasibility = current_iterate.progress.infeasibility;
      const bool h_type_iterates_allowed = this->constraint_relaxation_strategy->solving_feasibility_problem() ||
         this->constraint_relaxation_strategy->optimality_small_infeasibility_threshold() < infeasibility ||
         direction.norm <= this->h_type_step_bound(infeasibility);
      // restore the default (allowed) even if the acceptance test throws (e.g. evaluation error)
      struct HTypeIteratesFlagGuard {
         ConstraintRelaxationStrategy& strategy;
         ~HTypeIteratesFlagGuard() { this->strategy.set_h_type_iterates_allowed(true); }
      } guard{*this->constraint_relaxation_strategy};
      this->constraint_relaxation_strategy->set_h_type_iterates_allowed(h_type_iterates_allowed);
      return this->constraint_relaxation_strategy->is_iterate_acceptable(statistics, model, current_iterate,
         trial_iterate, direction, 1., true, evaluation_cache.current_evaluations, evaluation_cache.trial_evaluations,
         predicted_reductions, warmstart_information, user_callbacks);
   }

   // check whether the algorithm can terminate at an iterate (with computed residuals) after a rejected step with very
   // small norm
   bool TrustRegionStrategy::check_termination_with_small_step(Iterate& iterate) const {
      // terminate with a feasible point
      if (iterate.scaled_residuals.primal_infeasibility <= this->primal_tolerance) {
         iterate.status = SolutionStatus::FEASIBLE_SMALL_STEP;
         return true;
      }
      else if (this->constraint_relaxation_strategy->solving_feasibility_problem()) { // terminate with an infeasible point
         iterate.status = SolutionStatus::INFEASIBLE_SMALL_STEP;
         return true;
      }
      else { // do not terminate, infeasible non stationary
         return false;
      }
   }

   void TrustRegionStrategy::possibly_increase_radius(double step_norm) {
      // increase the radius if the trust-region is active
      if (step_norm >= this->radius - this->activity_tolerance) {
         this->radius *= this->increase_factor;
         DEBUG << "Trust-region radius increased to " << this->radius << '\n';
      }
   }

   void TrustRegionStrategy::decrease_radius(double step_norm) {
      // reduce the radius to a value smaller than the primal step norm (otherwise, the reduction won't have an effect)
      this->radius = std::min(this->radius, step_norm) / this->decrease_factor;
      DEBUG << "Trust-region radius decreased to " << this->radius << '\n';
   }

   void TrustRegionStrategy::decrease_radius() {
      this->radius /= this->decrease_factor;
      DEBUG << "Trust-region radius decreased to " << this->radius << '\n';
   }

   void TrustRegionStrategy::decrease_radius_aggressively() {
      this->radius /= this->aggressive_decrease_factor;
      DEBUG << "Trust-region radius aggressively decreased to " << this->radius << '\n';
   }

   void TrustRegionStrategy::reset_radius() {
      this->radius = std::max(this->radius, this->radius_reset_threshold);
   }

   void TrustRegionStrategy::set_TR_statistics(Statistics& statistics, size_t number_iterations) const {
      statistics.set("TR", number_iterations);
      statistics.set("Radius", this->radius);
   }
} // namespace