// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONSTRATEGY_H
#define UNO_GLOBALIZATIONSTRATEGY_H

namespace uno {
   // forward declarations
   class Iterate;
   struct ProgressMeasures;
   class Statistics;
   class Options;

   /*! \class GlobalizationStrategy
    *  Ingredient that accepts or rejects a trial iterate
    */
   class GlobalizationStrategy {
   public:
      explicit GlobalizationStrategy(const Options& options);
      virtual ~GlobalizationStrategy() = default;

      virtual void initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) = 0;
      [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
            const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) = 0;
      [[nodiscard]] virtual bool is_infeasibility_sufficiently_reduced(const ProgressMeasures& reference_progress, const ProgressMeasures& trial_progress) const = 0;

      virtual void reset() = 0;

      virtual void notify_switch_to_feasibility(const ProgressMeasures& current_progress) = 0;
      virtual void notify_switch_to_optimality(const ProgressMeasures& current_progress) = 0;

   protected:
      const double armijo_decrease_fraction; /*!< Sufficient reduction constant */
      const double armijo_tolerance;
      const bool protect_actual_reduction_against_roundoff;

      [[nodiscard]] bool armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const;
   };
} // namespace

#endif // UNO_GLOBALIZATIONSTRATEGY_H
