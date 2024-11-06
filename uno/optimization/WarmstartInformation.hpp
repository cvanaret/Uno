// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WARMSTARTINFORMATION_H
#define UNO_WARMSTARTINFORMATION_H

namespace uno {
   struct WarmstartInformation {
      bool objective_changed{false};
      bool constraints_changed{false};
      bool constraint_bounds_changed{false};
      bool variable_bounds_changed{false};
      bool problem_changed{false};

      void display() const;
      void set_cold_start();
      void set_hot_start();
      void only_objective_changed();
      void only_variable_bounds_changed();
   };
} // namespace

#endif // UNO_WARMSTARTINFORMATION_H
