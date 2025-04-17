// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "InequalityHandlingMethod.hpp"

namespace uno {
   void InequalityHandlingMethod::set_trust_region_radius(double new_trust_region_radius) {
      assert(0. < new_trust_region_radius && "The trust-region radius should be positive.");
      this->trust_region_radius = new_trust_region_radius;
   }
} // namespace