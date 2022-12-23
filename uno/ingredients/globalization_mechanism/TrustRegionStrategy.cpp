// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <cassert>
#include "TrustRegionStrategy.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"

TrustRegionStrategy::TrustRegionStrategy(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options) :
      GlobalizationMechanism(constraint_relaxation_strategy),
      radius(options.get_double("TR_radius")),
      increase_factor(options.get_double("TR_increase_factor")),
      decrease_factor(options.get_double("TR_decrease_factor")),
      activity_tolerance(options.get_double("TR_activity_tolerance")),
      min_radius(options.get_double("TR_min_radius")),
      statistics_TR_radius_column_order(options.get_int("statistics_TR_radius_column_order")) {
   assert(0 < this->radius && "The trust-region radius should be positive");
   assert(1. < this->increase_factor && "The trust-region increase factor should be > 1");
   assert(1. < this->decrease_factor && "The trust-region decrease factor should be > 1");
}

void TrustRegionStrategy::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("TR radius", Statistics::double_width, this->statistics_TR_radius_column_order);

   // generate the initial point
   this->constraint_relaxation_strategy.set_trust_region_radius(this->radius);
   this->constraint_relaxation_strategy.initialize(statistics, first_iterate);
}

std::tuple<Iterate, double> TrustRegionStrategy::compute_acceptable_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate) {
   this->number_iterations = 0;

   while (!this->termination()) {
      try {
         this->number_iterations++;
         this->print_iteration();

         // compute the direction within the trust region
         this->constraint_relaxation_strategy.set_trust_region_radius(this->radius);
         Direction direction = this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate);

         // assemble the trial iterate by taking a full step
         Iterate trial_iterate = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction);
         // reset bound multipliers of active trust region
         this->reset_active_trust_region_multipliers(model, direction, trial_iterate);

         // check whether the trial step is accepted
         const bool is_acceptable = this->constraint_relaxation_strategy
               .is_iterate_acceptable(statistics, current_iterate, trial_iterate, direction, 1.);
         if (is_acceptable) {
            // let the subproblem know the accepted iterate
            this->constraint_relaxation_strategy.register_accepted_iterate(trial_iterate);
            this->add_statistics(statistics, direction);

            // increase the radius if trust region is active
            this->increase_radius(direction.norm);

            return std::make_tuple(std::move(trial_iterate), direction.norm);
         }
         else { // step rejected
            this->decrease_radius(direction.norm);
         }
      }
      // if an error occurs (evaluation error or unstable inertia), decrease the radius
      catch (const std::exception& e) {
         GlobalizationMechanism::print_warning(e.what());
         this->decrease_radius();
      }
   }
   throw std::runtime_error("Trust-region radius became too small");
}

void TrustRegionStrategy::increase_radius(double step_norm) {
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

void TrustRegionStrategy::add_statistics(Statistics& statistics, const Direction& direction) {
   statistics.add_statistic("minor", this->number_iterations);
   statistics.add_statistic("TR radius", this->radius);
   statistics.add_statistic("step norm", direction.norm);
}

bool TrustRegionStrategy::termination() const {
   return this->radius < this->min_radius;
}

void TrustRegionStrategy::print_iteration() {
   DEBUG << "\t### Trust-region inner iteration " << this->number_iterations << " with radius " << this->radius << "\n\n";
}
