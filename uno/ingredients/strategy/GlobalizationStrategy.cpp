#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(const Options& options):
   armijo_decrease_fraction(stod(options.at("armijo_decrease_fraction"))),
   armijo_tolerance(stod(options.at("armijo_tolerance"))) {}

bool GlobalizationStrategy::armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const {
   return (actual_reduction >= this->armijo_decrease_fraction * std::max(0., predicted_reduction - this->armijo_tolerance));
}

void GlobalizationStrategy::check_finiteness(const ProgressMeasures& progress) {
   assert(is_finite(progress.optimality) && "The objective measure is infinite.");
   assert(is_finite(progress.infeasibility) && "The infeasibility measure is infinite.");
}