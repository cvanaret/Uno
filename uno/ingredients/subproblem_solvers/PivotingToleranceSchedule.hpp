// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PIVOTINGTOLERANCESCHEDULE_H
#define UNO_PIVOTINGTOLERANCESCHEDULE_H

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

namespace uno {
   // adaptive pivoting tolerance of a direct linear solver (as in IPOPT): u <- min(u_max, u^exponent).
   // Increases can be undone (in reverse order), down to the default tolerance
   class PivotingToleranceSchedule {
   public:
      PivotingToleranceSchedule(double maximum_tolerance, double increase_exponent):
            maximum_tolerance(maximum_tolerance), increase_exponent(increase_exponent) { }

      // returns the increased tolerance, or nothing if the maximum is reached
      [[nodiscard]] std::optional<double> increase(double current_tolerance) {
         if (this->maximum_tolerance <= current_tolerance) {
            return std::nullopt;
         }
         this->previous_tolerances.push_back(current_tolerance);
         return std::min(this->maximum_tolerance, std::pow(current_tolerance, this->increase_exponent));
      }

      // returns the tolerance before the last increase, or nothing if the tolerance was never increased
      [[nodiscard]] std::optional<double> decrease() {
         if (this->previous_tolerances.empty()) {
            return std::nullopt;
         }
         const double tolerance = this->previous_tolerances.back();
         this->previous_tolerances.pop_back();
         return tolerance;
      }

   protected:
      const double maximum_tolerance;
      const double increase_exponent;
      std::vector<double> previous_tolerances{};
   };
} // namespace

#endif // UNO_PIVOTINGTOLERANCESCHEDULE_H
