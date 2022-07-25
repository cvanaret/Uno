// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <cassert>
#include "TrustRegionStrategy.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"

TrustRegionStrategy::TrustRegionStrategy(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options) :
      GlobalizationMechanism(constraint_relaxation_strategy),
      radius(stod(options.at("TR_radius"))),
      increase_factor(stod(options.at("TR_increase_factor"))),
      decrease_factor(stod(options.at("TR_decrease_factor"))),
      activity_tolerance(stod(options.at("TR_activity_tolerance"))),
      min_radius(stod(options.at("TR_min_radius"))),
      statistics_TR_radius_column_order(stoi(options.at("statistics_TR_radius_column_order"))) {
   assert(1. < this->increase_factor && "The TR increase factor should be > 1");
   assert(1. < this->decrease_factor && "The TR decrease factor should be > 1");
}

void TrustRegionStrategy::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("TR radius", Statistics::double_width, this->statistics_TR_radius_column_order);

   // generate the initial point
   this->constraint_relaxation_strategy.set_variable_bounds(first_iterate, this->radius);
   this->constraint_relaxation_strategy.initialize(statistics, first_iterate);
}

std::tuple<Iterate, double> TrustRegionStrategy::compute_acceptable_iterate(Statistics& statistics, Iterate& current_iterate) {
   this->number_iterations = 0;

   while (!this->termination()) {
      try {
         assert (0 < this->radius);
         this->number_iterations++;
         this->print_iteration();

         // compute the direction within the trust region
         this->constraint_relaxation_strategy.set_variable_bounds(current_iterate, this->radius);
         Direction direction = this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate);
         // set bound multipliers of active trust region to 0
         TrustRegionStrategy::rectify_active_set(direction, this->radius);

         // assemble the trial iterate by taking a full step
         Iterate trial_iterate = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction);

         // check whether the trial step is accepted
         PredictedOptimalityReductionModel predicted_optimality_reduction_model = this->constraint_relaxation_strategy.generate_predicted_optimality_reduction_model(direction);
         if (this->constraint_relaxation_strategy.is_acceptable(statistics, current_iterate, trial_iterate, direction,
               predicted_optimality_reduction_model, 1.)) {
            DEBUG << "Trial step accepted\n\n";
            this->set_statistics(statistics, direction);

            // increase the radius if trust region is active
            if (direction.norm >= this->radius - this->activity_tolerance) {
               this->radius *= this->increase_factor;
            }

            // let the subproblem know the accepted iterate
            this->constraint_relaxation_strategy.register_accepted_iterate(trial_iterate);
            return std::make_tuple(std::move(trial_iterate), direction.norm);
         }
         else {
            // if the step is rejected, decrease the radius
            this->radius = std::min(this->radius, direction.norm) / this->decrease_factor;
         }
      }
      catch (const std::exception& e) {
         GlobalizationMechanism::print_warning(e.what());
         // if an error occurs (evaluation or unstable inertia), decrease the radius
         this->radius /= this->decrease_factor;
      }
   }
   // radius gets too small
   if (this->radius < this->min_radius) {
      throw std::runtime_error("Trust-region radius became too small");
   }
   else {
      throw std::runtime_error("Trust-region failed with an unexpected error");
   }
}

void TrustRegionStrategy::rectify_active_set(Direction& direction, double radius) {
   assert(0 < radius && "The trust-region radius is zero");
   // update active set and set multipliers for bound constraints active at trust region to 0
   for (auto it = direction.active_set.bounds.at_lower_bound.begin(); it != direction.active_set.bounds.at_lower_bound.end();) {
      size_t i = *it;
      if (direction.primals[i] == -radius) {
         it = direction.active_set.bounds.at_lower_bound.erase(it);
         direction.multipliers.lower_bounds[i] = 0.;
      }
      else {
         ++it;
      }
   }
   for (auto it = direction.active_set.bounds.at_upper_bound.begin(); it != direction.active_set.bounds.at_upper_bound.end();) {
      size_t i = *it;
      if (direction.primals[i] == radius) {
         it = direction.active_set.bounds.at_upper_bound.erase(it);
         direction.multipliers.upper_bounds[i] = 0.;
      }
      else {
         ++it;
      }
   }
}

void TrustRegionStrategy::set_statistics(Statistics& statistics, const Direction& direction) {
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
