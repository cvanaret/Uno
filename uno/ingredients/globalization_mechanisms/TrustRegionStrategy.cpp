// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <cassert>
#include "TrustRegionStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   TrustRegionStrategy::TrustRegionStrategy(const Options& options) :
         GlobalizationMechanism(),
         radius(options.get_double("TR_radius")),
         increase_factor(options.get_double("TR_increase_factor")),
         decrease_factor(options.get_double("TR_decrease_factor")),
         aggressive_decrease_factor(options.get_double("TR_aggressive_decrease_factor")),
         activity_tolerance(options.get_double("TR_activity_tolerance")),
         minimum_radius(options.get_double("TR_min_radius")),
         radius_reset_threshold(options.get_double("TR_radius_reset_threshold")),
         tolerance(options.get_double("tolerance")) {
      assert(0 < this->radius && "The trust-region radius should be positive");
      assert(1. < this->increase_factor && "The trust-region increase factor should be > 1");
      assert(1. < this->decrease_factor && "The trust-region decrease factor should be > 1");
   }

   void TrustRegionStrategy::initialize(Statistics& statistics, const Options& options) {
      statistics.add_column("TR iter", Statistics::int_width + 2, options.get_int("statistics_minor_column_order"));
      statistics.add_column("TR radius", Statistics::double_width - 4, options.get_int("statistics_TR_radius_column_order"));
      statistics.set("TR radius", this->radius);
   }

   void TrustRegionStrategy::compute_next_iterate(Statistics& statistics, ConstraintRelaxationStrategy& constraint_relaxation_strategy,
         GlobalizationStrategy& globalization_strategy, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         Direction& direction, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      DEBUG2 << "Current iterate\n" << current_iterate << '\n';
      this->reset_radius();

      size_t number_iterations = 0;
      bool termination = false;
      while (!termination) {
         bool is_acceptable = false;
         try {
            number_iterations++;
            DEBUG << "\n\t### Trust-region inner iteration " << number_iterations << " with radius " << this->radius << "\n\n";
            if (1 < number_iterations) { statistics.start_new_line(); }
            this->set_trust_region_statistics(statistics, number_iterations);

            // compute the direction within the trust region
            constraint_relaxation_strategy.compute_feasible_direction(statistics, globalization_strategy, model, current_iterate,
               direction, this->radius, warmstart_information);

            // deal with errors in the subproblem
            if (direction.status == SubproblemStatus::UNBOUNDED_PROBLEM) {
               // the subproblem is always bounded, but the objective may exceed a very large negative value
               this->set_statistics(statistics, direction);
               statistics.set("status", "unbounded subproblem");
               if (Logger::level == INFO) statistics.print_current_line();
               this->decrease_radius_aggressively();
               warmstart_information.variable_bounds_changed = true;
            }
            else if (direction.status == SubproblemStatus::ERROR) {
               this->set_statistics(statistics, direction);
               statistics.set("status", "solver error");
               if (Logger::level == INFO) statistics.print_current_line();
               this->decrease_radius();
               // reset the Hessian representation of the subproblem solver
               warmstart_information.whole_problem_changed();
            }
            else {
               // take full primal-dual step
               GlobalizationMechanism::assemble_trial_iterate(model, current_iterate, trial_iterate, direction, 1., 1.);
               this->reset_active_trust_region_multipliers(model, direction, trial_iterate);

               is_acceptable = this->is_iterate_acceptable(statistics, constraint_relaxation_strategy, globalization_strategy,
                  model, current_iterate, trial_iterate, direction, warmstart_information, user_callbacks);
               if (is_acceptable) {
                  constraint_relaxation_strategy.set_dual_residuals_statistics(statistics, trial_iterate);
                  termination = true;
               }
               else {
                  this->decrease_radius(direction.norm);
                  warmstart_information.variable_bounds_changed = true;
               }
               if (Logger::level == INFO) statistics.print_current_line();
            }
         }
         // if an evaluation error occurs, decrease the radius
         catch (const EvaluationError&) {
            statistics.set("status", "eval. error");
            if (Logger::level == INFO) statistics.print_current_line();
            DEBUG << "A function could not be evaluated. The trust-region radius will be reduced\n";
            this->decrease_radius();
            warmstart_information.variable_bounds_changed = true;
         }
         if (!is_acceptable && this->radius < this->minimum_radius) {
            throw std::runtime_error("Small radius");
         }
      }
   }

   std::string TrustRegionStrategy::get_name() const {
      return "TR";
   }

   // protected member functions

   void TrustRegionStrategy::reset_active_trust_region_multipliers(const Model& model, const Direction& direction, Iterate& trial_iterate) const {
      assert(0 < this->radius && "The trust-region radius should be positive");
      // reset multipliers for bound constraints active at trust region (except if one of the original bounds is active)
      for (size_t variable_index: Range(model.number_variables)) {
         if (std::abs(direction.primals[variable_index] + this->radius) <= this->activity_tolerance &&
               this->activity_tolerance < std::abs(trial_iterate.primals[variable_index] - model.variable_lower_bound(variable_index))) {
            trial_iterate.multipliers.lower_bounds[variable_index] = 0.;
            trial_iterate.feasibility_multipliers.lower_bounds[variable_index] = 0.;
         }
         if (std::abs(direction.primals[variable_index] - this->radius) <= this->activity_tolerance &&
               this->activity_tolerance < std::abs(model.variable_upper_bound(variable_index) - trial_iterate.primals[variable_index])) {
            trial_iterate.multipliers.upper_bounds[variable_index] = 0.;
            trial_iterate.feasibility_multipliers.upper_bounds[variable_index] = 0.;
         }
      }
   }

   // the trial iterate is accepted by the constraint relaxation strategy or if the step is small and we cannot switch to solving the feasibility problem
   bool TrustRegionStrategy::is_iterate_acceptable(Statistics& statistics, ConstraintRelaxationStrategy& constraint_relaxation_strategy,
         GlobalizationStrategy& globalization_strategy, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      bool accept_iterate = constraint_relaxation_strategy.is_iterate_acceptable(statistics, globalization_strategy, model,
         current_iterate, trial_iterate, direction, 1., warmstart_information, user_callbacks);
      this->set_statistics(statistics, trial_iterate, direction);
      if (accept_iterate) {
         trial_iterate.status = constraint_relaxation_strategy.check_termination(model, trial_iterate);
         // possibly increase the radius if trust region is active
         this->possibly_increase_radius(direction.norm);
      }
      else if (this->radius < this->minimum_radius) { // rejected, but small radius
         accept_iterate = this->check_termination_with_small_step(constraint_relaxation_strategy, model, trial_iterate);
      }
      return accept_iterate;
   }

   bool TrustRegionStrategy::check_termination_with_small_step(ConstraintRelaxationStrategy& constraint_relaxation_strategy,
         const Model& model, Iterate& trial_iterate) const {
      // terminate with a feasible point
      if (trial_iterate.progress.infeasibility <= this->tolerance) {
         trial_iterate.status = IterateStatus::FEASIBLE_SMALL_STEP;
         constraint_relaxation_strategy.compute_primal_dual_residuals(model, trial_iterate);
         return true;
      }
      else if (constraint_relaxation_strategy.solving_feasibility_problem()) { // terminate with an infeasible point
         trial_iterate.status = IterateStatus::INFEASIBLE_SMALL_STEP;
         constraint_relaxation_strategy.compute_primal_dual_residuals(model, trial_iterate);
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

   void TrustRegionStrategy::set_trust_region_statistics(Statistics& statistics, size_t number_iterations) const {
      statistics.set("TR iter", number_iterations);
      statistics.set("TR radius", this->radius);
   }

   void TrustRegionStrategy::set_statistics(Statistics& statistics, const Direction& direction) const {
      statistics.set("step norm", direction.norm);
   }

   void TrustRegionStrategy::set_statistics(Statistics& statistics, const Iterate& trial_iterate, const Direction& direction) const {
      if (trial_iterate.is_objective_computed) {
         statistics.set("objective", trial_iterate.evaluations.objective);
      }
      this->set_statistics(statistics, direction);
   }
} // namespace
