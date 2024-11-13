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
      void no_changes();
      void iterate_changed();
      void whole_problem_changed();
      void only_objective_changed();
   };
} // namespace

#endif // UNO_WARMSTARTINFORMATION_H
