// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <cassert>
#include "TrustRegionStrategy.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"

TrustRegionStrategy::TrustRegionStrategy(Statistics& statistics, ConstraintRelaxationStrategy& constraint_relaxation_strategy,
         const Options& options) :
      GlobalizationMechanism(constraint_relaxation_strategy, options),
      radius(options.get_double("TR_radius")),
      increase_factor(options.get_double("TR_increase_factor")),
      decrease_factor(options.get_double("TR_decrease_factor")),
      aggressive_decrease_factor(options.get_double("TR_aggressive_decrease_factor")),
      activity_tolerance(options.get_double("TR_activity_tolerance")),
      minimum_radius(options.get_double("TR_min_radius")),
      radius_reset_threshold(options.get_double("TR_radius_reset_threshold")) {
   assert(0 < this->radius && "The trust-region radius should be positive");
   assert(1. < this->increase_factor && "The trust-region increase factor should be > 1");
   assert(1. < this->decrease_factor && "The trust-region decrease factor should be > 1");

   statistics.add_column("TR iters", Statistics::int_width + 3, options.get_int("statistics_minor_column_order"));
   statistics.add_column("TR radius", Statistics::double_width, options.get_int("statistics_TR_radius_column_order"));
}

void TrustRegionStrategy::initialize(Iterate& initial_iterate) {
   this->constraint_relaxation_strategy.set_trust_region_radius(this->radius);
   this->constraint_relaxation_strategy.initialize(initial_iterate);
}

Iterate TrustRegionStrategy::compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate) {
   this->number_iterations = 0;
   this->reset_radius();

   WarmstartInformation warmstart_information{};
   warmstart_information.objective_changed = true;
   warmstart_information.constraints_changed = true;
   warmstart_information.constraint_bounds_changed = true;
   warmstart_information.variable_bounds_changed = true;
   DEBUG2 << "Current iterate\n" << current_iterate << '\n';

   while (not this->termination()) {
      try {
         this->number_iterations++;
         this->print_iteration();

         // compute the direction within the trust region
         this->constraint_relaxation_strategy.set_trust_region_radius(this->radius);
         Direction direction = this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate, warmstart_information);
         DEBUG << "Step norm: " << direction.norm << '\n';

         if (direction.status == SubproblemStatus::UNBOUNDED_PROBLEM) {
            this->decrease_radius_aggressively();
            warmstart_information.objective_changed = true;
            warmstart_information.constraints_changed = true;
            warmstart_information.constraint_bounds_changed = true;
            warmstart_information.variable_bounds_changed = true;
            warmstart_information.problem_changed = true;
         }
         if (direction.status == SubproblemStatus::ERROR) {
            this->decrease_radius();
            warmstart_information.objective_changed = true;
            warmstart_information.constraints_changed = true;
            warmstart_information.constraint_bounds_changed = true;
            warmstart_information.variable_bounds_changed = true;
            warmstart_information.problem_changed = true;
         }
         else {
            // assemble the trial iterate by taking a full step
            Iterate trial_iterate = this->assemble_trial_iterate(model, current_iterate, direction);
            // check whether the trial step is accepted
            if (this->constraint_relaxation_strategy.is_iterate_acceptable(statistics, current_iterate, trial_iterate, direction,
                  direction.primal_dual_step_length)) {
               this->set_statistics(statistics, direction);

               // increase the radius if trust region is active
               this->possibly_increase_radius(direction.norm);

               // check termination criteria
               trial_iterate.status = this->check_termination(model, trial_iterate, direction.norm);
               return trial_iterate;
            }
            else { // trial iterate not acceptable
               this->decrease_radius(direction.norm);
            }
            // after the first iteration, only the variable bounds are updated
            warmstart_information.objective_changed = false;
            warmstart_information.constraints_changed = false;
            warmstart_information.constraint_bounds_changed = false;
            warmstart_information.variable_bounds_changed = true;
         }
      }
      catch (const std::runtime_error& e) {
         throw;
      }
      // if an error occurs (evaluation error or unstable inertia), decrease the radius
      catch (const std::exception& e) {
         WARNING << YELLOW << e.what() << RESET;
         this->decrease_radius();
         warmstart_information.objective_changed = true;
         warmstart_information.constraints_changed = true;
         warmstart_information.constraint_bounds_changed = true;
         warmstart_information.variable_bounds_changed = true;
      }
   }
   // TODO: may still be accepted as solution if the termination criteria are satisfied
   throw std::runtime_error("Trust-region radius became too small\n");
}

Iterate TrustRegionStrategy::assemble_trial_iterate(const Model& model, Iterate& current_iterate, const Direction& direction) {
   Iterate trial_iterate = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction, direction.primal_dual_step_length,
         direction.bound_dual_step_length);
   // project the steps within the bounds to avoid numerical errors
   model.project_primals_onto_bounds(trial_iterate.primals);

   // reset bound multipliers of active trust region
   this->reset_active_trust_region_multipliers(model, direction, trial_iterate);
   return trial_iterate;
}

void TrustRegionStrategy::possibly_increase_radius(double step_norm) {
   if (step_norm >= this->radius - this->activity_tolerance) {
      this->radius *= this->increase_factor;
   }
}

void TrustRegionStrategy::decrease_radius(double step_norm) {
   this->radius = std::min(this->radius, step_norm) / this->decrease_factor;
}

void TrustRegionStrategy::decrease_radius() {
   this->radius /= this->decrease_factor;
}

void TrustRegionStrategy::decrease_radius_aggressively() {
   this->radius /= this->aggressive_decrease_factor;
}

void TrustRegionStrategy::reset_radius() {
   this->radius = std::max(this->radius, this->radius_reset_threshold);
}

void TrustRegionStrategy::reset_active_trust_region_multipliers(const Model& model, const Direction& direction, Iterate& trial_iterate) const {
   assert(0 < this->radius && "The trust-region radius should be positive");
   // set multipliers for bound constraints active at trust region to 0 (except if one of the original bounds is active)
   for (size_t i: direction.active_set.bounds.at_lower_bound) {
      if (i < model.number_variables && std::abs(direction.primals[i] + this->radius) <= this->activity_tolerance &&
            this->activity_tolerance < std::abs(trial_iterate.primals[i] - model.get_variable_lower_bound(i))) {
         trial_iterate.multipliers.lower_bounds[i] = 0.;
      }
   }
   for (size_t i: direction.active_set.bounds.at_upper_bound) {
      if (i < model.number_variables && std::abs(direction.primals[i] - this->radius) <= this->activity_tolerance &&
            this->activity_tolerance < std::abs(model.get_variable_upper_bound(i) - trial_iterate.primals[i])) {
         trial_iterate.multipliers.upper_bounds[i] = 0.;
      }
   }
}

void TrustRegionStrategy::set_statistics(Statistics& statistics, const Direction& direction) {
   statistics.add_statistic("TR iters", this->number_iterations);
   statistics.add_statistic("TR radius", this->radius);
   statistics.add_statistic("step norm", direction.norm);
}

bool TrustRegionStrategy::termination() const {
   return this->radius < this->minimum_radius;
}

void TrustRegionStrategy::print_iteration() {
   DEBUG << "\t### Trust-region inner iteration " << this->number_iterations << " with radius " << this->radius << "\n\n";
}
