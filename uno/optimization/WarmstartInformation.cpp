// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include "WarmstartInformation.hpp"

namespace uno {
   void WarmstartInformation::no_changes() {
      this->iterate_changed = false;
      this->variable_bounds_changed = false;
   }

   void WarmstartInformation::reset() {
      this->iterate_changed = true;
      this->variable_bounds_changed = true;
   }

   void WarmstartInformation::display() const {
      std::cout << "Iterate changed: " << std::boolalpha << this->iterate_changed << '\n';
      std::cout << "Variable bounds changed: " << std::boolalpha << this->variable_bounds_changed << '\n';
   }
} // namespace