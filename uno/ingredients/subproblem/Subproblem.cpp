// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "Subproblem.hpp"

Subproblem::Subproblem(size_t max_number_variables, size_t max_number_constraints):
      direction(max_number_variables, max_number_constraints),
      evaluations(max_number_variables, max_number_constraints) {
}

void Subproblem::set_trust_region_radius(double new_trust_region_radius) {
   assert(0. < new_trust_region_radius && "The trust-region radius should be positive.");
   this->trust_region_radius = new_trust_region_radius;
}

void Subproblem::check_unboundedness([[maybe_unused]] const Direction& direction) {
   assert(direction.status != SubproblemStatus::UNBOUNDED_PROBLEM && "The subproblem is unbounded, this should not happen");
}
