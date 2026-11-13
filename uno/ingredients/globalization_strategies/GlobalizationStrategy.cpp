// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
#include <limits>
#include "GlobalizationStrategy.hpp"
#include "options/Options.hpp"

namespace uno {
   GlobalizationStrategy::GlobalizationStrategy(const Options& options):
      armijo_decrease_fraction(options.get_double("armijo_decrease_fraction")),
      armijo_tolerance(options.get_double("armijo_tolerance")),
      protect_actual_reduction_against_roundoff(options.get_bool("protect_actual_reduction_against_roundoff")),
      protected_actual_reduction_macheps_coefficient(options.get_double("protected_actual_reduction_macheps_coefficient")),
      small_infeasibility_factor(options.get_double("barrier_small_infeasibility_factor")) {
   }

   GlobalizationStrategy::~GlobalizationStrategy() = default;

   // default: absolute threshold
   double GlobalizationStrategy::small_infeasibility_threshold() const {
      return this->small_infeasibility_factor;
   }

   bool GlobalizationStrategy::armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const {
      return (actual_reduction >= this->armijo_decrease_fraction * std::max(0., predicted_reduction - this->armijo_tolerance));
   }

   // actual reduction of a progress measure, possibly protected against roundoff errors
   double GlobalizationStrategy::compute_actual_reduction(double current_measure, double trial_measure) const {
      double actual_reduction = current_measure - trial_measure;
      if (this->protect_actual_reduction_against_roundoff) {
         constexpr double machine_epsilon = std::numeric_limits<double>::epsilon();
         actual_reduction += this->protected_actual_reduction_macheps_coefficient * machine_epsilon * std::abs(current_measure);
      }
      return actual_reduction;
   }
} // namespace