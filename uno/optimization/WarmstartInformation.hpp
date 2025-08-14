// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WARMSTARTINFORMATION_H
#define UNO_WARMSTARTINFORMATION_H

namespace uno {
   class WarmstartInformation {
   public:
      bool objective_changed{true};
      bool constraints_changed{true};
      bool constraint_bounds_changed{true};
      bool variable_bounds_changed{true};

      void no_changes();
      void reset();

      void display() const;
   };
} // namespace

#endif // UNO_WARMSTARTINFORMATION_H
