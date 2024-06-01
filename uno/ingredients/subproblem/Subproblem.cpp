// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "Subproblem.hpp"

Subproblem::Subproblem(size_t number_variables, size_t number_constraints):
      evaluations(number_variables, number_constraints) {
}

void Subproblem::set_trust_region_radius(double new_trust_region_radius) {
   assert(0. < new_trust_region_radius && "The trust-region radius should be positive.");
   this->trust_region_radius = new_trust_region_radius;
}