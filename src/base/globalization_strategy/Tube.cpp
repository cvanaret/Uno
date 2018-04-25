#include <iostream>
#include "Tube.hpp"
#include "Utils.hpp"

Tube::Tube(double upper_bound, StepConstants& constants): upper_bound(upper_bound), constants(constants) {
}

//! update upper bound on tube
void Tube::update(double current_residual, double trial_residual) {
	double actual_reduction = max(0.001, current_residual - trial_residual);
	this->upper_bound = max(0.9*this->upper_bound, this->upper_bound - 0.1*actual_reduction);
}

//! check whether new constraint violation is acceptable
bool Tube::query(double current_residual) {
	// acceptablity depends only on upper bound
	return (current_residual < this->upper_bound);
}
