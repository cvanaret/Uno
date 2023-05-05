// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WARMSTARTINFORMATION_H
#define UNO_WARMSTARTINFORMATION_H

struct WarmstartInformation {
   bool objective_changed{false};
   bool constraints_changed{false};
   bool constraint_bounds_changed{false};
   bool variable_bounds_changed{false};

   void display() const;
};

inline void WarmstartInformation::display() const {
   std::cout << "Objective: " << std::boolalpha << this->objective_changed << '\n';
   std::cout << "Constraints: " << std::boolalpha << this->constraints_changed << '\n';
   std::cout << "Constraint bounds: " << std::boolalpha << this->constraint_bounds_changed << '\n';
   std::cout << "Variable bounds: " << std::boolalpha << this->variable_bounds_changed << '\n';
}

#endif // UNO_WARMSTARTINFORMATION_H
