// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_REJECTIONHISTORY_H
#define UNO_REJECTIONHISTORY_H

#include "ingredients/globalization_strategies/Decision.hpp"

namespace uno {
   // rejections of the trial iterates of a single step (e.g. a backtracking line search). Trial iterates whose
   // evaluation failed were not tested, and are not recorded
   class RejectionHistory {
   public:
      void record(const Decision& decision) {
         if (decision.accepted) {
            return;
         }
         if (this->first_reason == RejectionReason::NONE) {
            this->first_reason = decision.reason;
         }
         else if (decision.reason != this->first_reason) {
            this->single_reason = false;
         }
      }

      // the reason shared by all the rejected trial iterates. NONE if no trial iterate was rejected or if the reasons differ
      [[nodiscard]] RejectionReason sole_reason() const {
         return this->single_reason ? this->first_reason : RejectionReason::NONE;
      }

   private:
      RejectionReason first_reason{RejectionReason::NONE};
      bool single_reason{true};
   };
} // namespace

#endif // UNO_REJECTIONHISTORY_H
