// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include "WarmstartInformation.hpp"

namespace uno {
   void WarmstartInformation::display() const {
      std::cout << "Objective: " << std::boolalpha << this->objective_changed << '\n';
      std::cout << "Constraints: " << std::boolalpha << this->constraints_changed << '\n';
      std::cout << "Constraint bounds: " << std::boolalpha << this->constraint_bounds_changed << '\n';
      std::cout << "Variable bounds: " << std::boolalpha << this->variable_bounds_changed << '\n';
      std::cout << "Problem: " << std::boolalpha << this->problem_changed << '\n';
   }

   void WarmstartInformation::no_changes() {
      this->objective_changed = false;
      this->constraints_changed = false;
      this->constraint_bounds_changed = false;
      this->variable_bounds_changed = false;
      this->problem_changed = false;
   }

   void WarmstartInformation::iterate_changed() {
      this->objective_changed = true;
      this->constraints_changed = true;
      this->constraint_bounds_changed = true;
      this->variable_bounds_changed = true;
   }

   void WarmstartInformation::whole_problem_changed() {
      this->objective_changed = true;
      this->constraints_changed = true;
      this->constraint_bounds_changed = true;
      this->variable_bounds_changed = true;
      this->problem_changed = true;
   }

   void WarmstartInformation::only_objective_changed() {
      this->objective_changed = true;
      this->constraints_changed = false;
      this->constraint_bounds_changed = false;
      this->variable_bounds_changed = false;
      this->problem_changed = false;
   }
} // namespace