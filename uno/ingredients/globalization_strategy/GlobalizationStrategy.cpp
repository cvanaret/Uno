// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(const Options& options):
   armijo_decrease_fraction(options.get_double("armijo_decrease_fraction")),
   armijo_tolerance(options.get_double("armijo_tolerance")) {}

bool GlobalizationStrategy::armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const {
   return (actual_reduction >= this->armijo_decrease_fraction * std::max(0., predicted_reduction - this->armijo_tolerance));
}

void GlobalizationStrategy::check_finiteness(const ProgressMeasures& progress) {
   assert(is_finite(progress.infeasibility) && "The infeasibility measure is infinite.");
   assert(is_finite(progress.scaled_optimality) && "The optimality measure is infinite.");
}