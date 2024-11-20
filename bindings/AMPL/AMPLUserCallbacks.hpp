// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_AMPLUSERCALLBACKS_H
#define UNO_AMPLUSERCALLBACKS_H

#include "tools/UserCallbacks.hpp"

namespace uno {
   class AMPLUserCallbacks: public UserCallbacks {
   public:
      AMPLUserCallbacks();

      void notify_acceptable_iterate(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier) override;
      void notify_new_primals(const Vector<double>& primals) override;
      void notify_new_multipliers(const Multipliers& multipliers) override;
   };
} // namespace

#endif //UNO_AMPLUSERCALLBACKS_H