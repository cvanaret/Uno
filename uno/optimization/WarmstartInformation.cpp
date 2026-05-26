// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include "WarmstartInformation.hpp"

namespace uno {
   void WarmstartInformation::display() const {
      std::cout << "New iterate: " << std::boolalpha << this->new_iterate << '\n';
      std::cout << "Constraint bounds changed: " << std::boolalpha << this->constraint_bounds_changed << '\n';
      std::cout << "Trust-region radius changed: " << std::boolalpha << this->trust_region_changed << '\n';
      std::cout << "Hessian sparsity changed: " << std::boolalpha << this->hessian_sparsity_changed << '\n';
      std::cout << "Jacobian sparsity changed: " << std::boolalpha << this->jacobian_sparsity_changed << '\n';
   }

   void WarmstartInformation::no_changes() {
      this->new_iterate = false;
      this->constraint_bounds_changed = false;
      this->trust_region_changed = false;
      this->hessian_sparsity_changed = false;
      this->jacobian_sparsity_changed = false;
   }

   void WarmstartInformation::iterate_changed() {
      this->new_iterate = true;
      this->constraint_bounds_changed = true;
      this->trust_region_changed = true;
      this->hessian_sparsity_changed = false;
      this->jacobian_sparsity_changed = false;
   }

   void WarmstartInformation::whole_problem_changed() {
      this->new_iterate = true;
      this->constraint_bounds_changed = true;
      this->trust_region_changed = true;
      this->hessian_sparsity_changed = true;
      this->jacobian_sparsity_changed = true;
   }

   void WarmstartInformation::only_objective_changed() {
      this->new_iterate = true;
      this->constraint_bounds_changed = false;
      this->trust_region_changed = false;
      this->hessian_sparsity_changed = false;
      this->jacobian_sparsity_changed = false;
   }
} // namespace