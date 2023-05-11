// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WARMSTARTINFORMATION_H
#define UNO_WARMSTARTINFORMATION_H

struct WarmstartInformation {
   bool objective_changed{false};
   bool constraints_changed{false};
   bool constraint_bounds_changed{false};
   bool variable_bounds_changed{false};
   bool problem_changed{false};

   void display() const;
   void set_cold_start();
   void set_hot_start();
   void only_variable_bounds_changed();
};

inline void WarmstartInformation::display() const {
   std::cout << "Objective: " << std::boolalpha << this->objective_changed << '\n';
   std::cout << "Constraints: " << std::boolalpha << this->constraints_changed << '\n';
   std::cout << "Constraint bounds: " << std::boolalpha << this->constraint_bounds_changed << '\n';
   std::cout << "Variable bounds: " << std::boolalpha << this->variable_bounds_changed << '\n';
   std::cout << "Problem: " << std::boolalpha << this->problem_changed << '\n';
}

inline void WarmstartInformation::set_cold_start() {
   this->objective_changed = true;
   this->constraints_changed = true;
   this->constraint_bounds_changed = true;
   this->variable_bounds_changed = true;
   this->problem_changed = true;
}

inline void WarmstartInformation::set_hot_start() {
   this->objective_changed = true;
   this->constraints_changed = true;
   this->constraint_bounds_changed = true;
   this->variable_bounds_changed = true;
   this->problem_changed = false;
}

inline void WarmstartInformation::only_variable_bounds_changed() {
   this->objective_changed = false;
   this->constraints_changed = false;
   this->constraint_bounds_changed = false;
   this->variable_bounds_changed = true;
   this->problem_changed = false;
}

#endif // UNO_WARMSTARTINFORMATION_H
