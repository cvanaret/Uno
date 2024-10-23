// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "GlobalizationStrategy.hpp"
#include "options/Options.hpp"

namespace uno {
   GlobalizationStrategy::GlobalizationStrategy(const Options& options):
      armijo_decrease_fraction(options.get_double("armijo_decrease_fraction")),
      armijo_tolerance(options.get_double("armijo_tolerance")),
      protect_actual_reduction_against_roundoff(options.get_bool("protect_actual_reduction_against_roundoff")) {
   }

   bool GlobalizationStrategy::armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const {
      return (actual_reduction >= this->armijo_decrease_fraction * std::max(0., predicted_reduction - this->armijo_tolerance));
   }
} // namespace
