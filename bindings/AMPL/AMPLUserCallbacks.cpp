// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "AMPLUserCallbacks.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Multipliers.hpp"

namespace uno {
   AMPLUserCallbacks::AMPLUserCallbacks(): UserCallbacks() { }

   void AMPLUserCallbacks::notify_acceptable_iterate(const Vector<double>& /*primals*/, const Multipliers& /*multipliers*/,
         double /*objective_multiplier*/) {
   }

   void AMPLUserCallbacks::notify_new_primals(const Vector<double>& /*primals*/) {
   }

   void AMPLUserCallbacks::notify_new_multipliers(const Multipliers& /*multipliers*/) {
   }
} // namespace