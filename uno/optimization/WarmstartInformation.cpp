// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include "WarmstartInformation.hpp"

namespace uno {
   void WarmstartInformation::no_changes() {
      this->objective_changed = false;
      this->constraints_changed = false;
      this->constraint_bounds_changed = false;
      this->variable_bounds_changed = false;
   }

   void WarmstartInformation::reset() {
      this->objective_changed = true;
      this->constraints_changed = true;
      this->constraint_bounds_changed = true;
      this->variable_bounds_changed = true;
   }

   void WarmstartInformation::display() const {
      std::cout << "Objective changed: " << std::boolalpha << this->objective_changed << '\n';
      std::cout << "Constraints changed: " << std::boolalpha << this->constraints_changed << '\n';
      std::cout << "Constraint bounds changed: " << std::boolalpha << this->constraint_bounds_changed << '\n';
      std::cout << "Variable bounds changed: " << std::boolalpha << this->variable_bounds_changed << '\n';
   }
} // namespace