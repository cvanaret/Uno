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

void GlobalizationStrategy::check_finiteness([[maybe_unused]] const ProgressMeasures& progress, [[maybe_unused]] double objective_multiplier) {
   assert(not std::isnan(progress.infeasibility) && is_finite(progress.infeasibility) && "The infeasibility measure is not finite.");
   assert(not std::isnan(progress.optimality(objective_multiplier)) && is_finite(progress.optimality(objective_multiplier)) && "The optimality measure is not finite.");
   assert(not std::isnan(progress.auxiliary_terms) && "The auxiliary measure is not a number.");
}
