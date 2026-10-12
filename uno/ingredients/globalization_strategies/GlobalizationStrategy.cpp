// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
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
} // namespace