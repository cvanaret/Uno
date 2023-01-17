// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(const Options& options):
   armijo_decrease_fraction(options.get_double("armijo_decrease_fraction")),
   armijo_tolerance(options.get_double("armijo_tolerance")) {}

bool GlobalizationStrategy::armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const {
   return (actual_reduction >= this->armijo_decrease_fraction * std::max(0., predicted_reduction - this->armijo_tolerance));
}

void GlobalizationStrategy::check_finiteness(const ProgressMeasures& progress, double objective_multiplier) {
   assert(not std::isnan(progress.infeasibility) && is_finite(progress.infeasibility) && "The infeasibility measure is not finite.");
   assert(not std::isnan(progress.scaled_optimality(objective_multiplier)) && is_finite(progress.scaled_optimality(objective_multiplier)) && "The scaled optimality measure is not finite.");
   assert(not std::isnan(progress.unscaled_optimality) && "The unscaled optimality measure is not a number.");
   //assert(is_finite(progress.unscaled_optimality) && "The unscaled optimality measure is infinite.");
}
