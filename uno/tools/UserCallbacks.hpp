// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_USERCALLBACKS_H
#define UNO_USERCALLBACKS_H

namespace uno {
   // forward declarations
   class Multipliers;
   template <class ElementType>
   class Vector;

   class UserCallbacks {
   public:
      UserCallbacks() = default;
      virtual ~UserCallbacks() = default;

      virtual void notify_acceptable_iterate(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier,
         double primal_feasibility_residual, double stationarity_residual, double complementarity_residual) = 0;
      virtual bool user_termination(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier,
         double primal_feasibility_residual, double stationarity_residual, double complementarity_residual) = 0; // returns true for user termination
   };

   class NoUserCallbacks: public UserCallbacks {
   public:
      NoUserCallbacks(): UserCallbacks() { }

      void notify_acceptable_iterate(const Vector<double>& /*primals*/, const Multipliers& /*multipliers*/, double /*objective_multiplier*/,
         double /*primal_feasibility_residual*/, double /*stationarity_residual*/, double /*complementarity_residual*/) override { }
      bool user_termination(const Vector<double>& /*primals*/, const Multipliers& /*multipliers*/, double /*objective_multiplier*/,
            double /*primal_feasibility_residual*/, double /*stationarity_residual*/, double /*complementarity_residual*/) override {
         return false;
      }
   };
} // namespace

#endif //UNO_USERCALLBACKS_H