#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(const Options& options):
   armijo_decrease_fraction(stod(options.at("armijo_decrease_fraction"))),
   armijo_tolerance(stod(options.at("armijo_tolerance"))) {}

bool GlobalizationStrategy::armijo_condition(double predicted_reduction, double actual_reduction) const {
   return actual_reduction >= this->armijo_decrease_fraction * std::max(0., predicted_reduction - this->armijo_tolerance);
}