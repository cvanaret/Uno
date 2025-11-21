// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONMECHANISM_H
#define UNO_GLOBALIZATIONMECHANISM_H

#include <memory>
#include <string>
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class Iterate;
   class Model;
   class Options;
   class Statistics;
   class UserCallbacks;
   class WarmstartInformation;

   class GlobalizationMechanism {
   public:
      GlobalizationMechanism(const Model& model, bool use_trust_region, const Options& options);
      virtual ~GlobalizationMechanism() = default;

      virtual void initialize(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Direction& direction) = 0;
      virtual void compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         Direction& direction, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) = 0;

      static void set_primal_statistics(Statistics& statistics, const Model& model, const Iterate& iterate);
      static void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate);

      [[nodiscard]] virtual std::string get_name() const = 0;
      [[nodiscard]] size_t get_number_subproblems_solved() const;

   protected:
      const std::unique_ptr<ConstraintRelaxationStrategy> constraint_relaxation_strategy{};

      static void assemble_trial_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double primal_step_length, double dual_step_length);
   };
} // namespace

#endif // UNO_GLOBALIZATIONMECHANISM_H