// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONMECHANISM_H
#define UNO_GLOBALIZATIONMECHANISM_H

#include <string>

namespace uno {
   // forward declarations
   class ConstraintRelaxationStrategy;
   class Direction;
   class GlobalizationStrategy;
   class Iterate;
   class Model;
   class Options;
   class Statistics;
   class UserCallbacks;
   class WarmstartInformation;

   class GlobalizationMechanism {
   public:
      GlobalizationMechanism() = default;
      virtual ~GlobalizationMechanism() = default;

      virtual void initialize(Statistics& statistics, const Options& options) = 0;
      virtual void compute_next_iterate(Statistics& statistics, ConstraintRelaxationStrategy& constraint_relaxation_strategy,
         GlobalizationStrategy& globalization_strategy, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         Direction& direction, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;

   protected:
      static void assemble_trial_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double primal_step_length, double dual_step_length);
   };
} // namespace

#endif // UNO_GLOBALIZATIONMECHANISM_H