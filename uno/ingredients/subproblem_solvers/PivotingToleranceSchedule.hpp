// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PIVOTINGTOLERANCESCHEDULE_H
#define UNO_PIVOTINGTOLERANCESCHEDULE_H

#include <algorithm>
#include <cmath>
#include <optional>

namespace uno {
   // adaptive pivoting tolerance of a direct linear solver (as in IPOPT): u <- min(u_max, u^exponent). The tolerance is
   // only ever increased (an increase is never undone)
   class PivotingToleranceSchedule {
   public:
      PivotingToleranceSchedule(double maximum_tolerance, double increase_exponent):
            maximum_tolerance(maximum_tolerance), increase_exponent(increase_exponent) { }

      // returns the increased tolerance, or nothing if the maximum is reached
      [[nodiscard]] std::optional<double> increase(double current_tolerance) const {
         if (this->maximum_tolerance <= current_tolerance) {
            return std::nullopt;
         }
         return std::min(this->maximum_tolerance, std::pow(current_tolerance, this->increase_exponent));
      }

   protected:
      const double maximum_tolerance;
      const double increase_exponent;
   };
} // namespace

#endif // UNO_PIVOTINGTOLERANCESCHEDULE_H
