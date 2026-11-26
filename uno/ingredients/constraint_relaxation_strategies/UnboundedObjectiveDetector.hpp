// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNBOUNDEDOBJECTIVEDETECTOR_H
#define UNO_UNBOUNDEDOBJECTIVEDETECTOR_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>

namespace uno {
   // Termination heuristic: the problem is declared unbounded when feasible iterates (primal infeasibility below the
   // tight unscaled tolerance, tested by the caller) have an objective that is
   // - below an absolute threshold (unbounded_objective_threshold), or
   // - lower than the objective f_ref of the first feasible iterate by many orders of magnitude:
   //   f < f_ref - relative_decrease * max(1, |f_ref|), at iteration_threshold consecutive feasible iterates.
   // On a feasible ray along which the objective decreases polynomially (CUTE coshfun: f = 5 - 2t^2), the absolute
   // threshold is out of reach within any iteration budget. No finite threshold certifies unboundedness: this only
   // decides when to stop, and leaves the iterates untouched
   class UnboundedObjectiveDetector {
   public:
      UnboundedObjectiveDetector(double absolute_threshold, double relative_decrease, size_t iteration_threshold):
            absolute_threshold(absolute_threshold),
            relative_decrease(relative_decrease),
            iteration_threshold(std::max(iteration_threshold, size_t(1))) {
      }

      // unscaled objective of a feasible iterate. Returns true if the problem is declared unbounded
      [[nodiscard]] bool is_unbounded(double objective) {
         if (objective < this->absolute_threshold) {
            return true;
         }
         if (!this->reference_objective.has_value()) {
            this->reference_objective = objective;
            return false;
         }
         const double reference = *this->reference_objective;
         const double threshold = reference - this->relative_decrease * std::max(1., std::abs(reference));
         if (objective < threshold) {
            ++this->streak;
            return (this->iteration_threshold <= this->streak);
         }
         this->streak = 0;
         return false;
      }

      // the streak counts consecutive feasible iterates
      void notify_infeasible_iterate() {
         this->streak = 0;
      }

   private:
      const double absolute_threshold;
      const double relative_decrease; // infinite: relative test disabled
      const size_t iteration_threshold;
      std::optional<double> reference_objective{};
      size_t streak{0};
   };
} // namespace

#endif // UNO_UNBOUNDEDOBJECTIVEDETECTOR_H
