// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DECISION_H
#define UNO_DECISION_H

namespace uno {
   // why a trial iterate was rejected
   enum class RejectionReason {
      NONE, // accepted
      INFEASIBILITY_UPPER_BOUND,
      SUFFICIENT_DECREASE, // Armijo condition
      CURRENT_ITERATE, // not acceptable wrt the current iterate
      FILTER,
      FUNNEL,
      INVALID_TRIAL_ITERATE // e.g. infinite progress measure
   };

   // acceptance decision on a trial iterate
   struct Decision {
      bool accepted{false};
      RejectionReason reason{RejectionReason::NONE};

      [[nodiscard]] static Decision accept() {
         return {true, RejectionReason::NONE};
      }
      [[nodiscard]] static Decision reject(RejectionReason reason) {
         return {false, reason};
      }
   };
} // namespace

#endif // UNO_DECISION_H
