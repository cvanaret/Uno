// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONSTRATEGY_H
#define UNO_GLOBALIZATIONSTRATEGY_H

#include <string>

namespace uno {
   // forward declarations
   class Iterate;
   class ProgressMeasures;
   class Statistics;
   class Options;

   class GlobalizationStrategy {
   public:
      explicit GlobalizationStrategy(const Options& options);
      virtual ~GlobalizationStrategy();

      virtual void initialize(Statistics& statistics, const Iterate& initial_iterate) = 0;
      [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reductions, double objective_multiplier) = 0;
      [[nodiscard]] virtual bool is_infeasibility_sufficiently_reduced(const Iterate& trial_iterate,
         double reference_infeasibility) const = 0;

      virtual void reset() = 0;

      // theta_min: below this infeasibility, the switching rule may declare f-type iterations. The line search uses the
      // same threshold to compute its minimum step length (IPOPT, Wächter & Biegler 2006, eq. (23))
      [[nodiscard]] virtual double small_infeasibility_threshold() const;

      virtual void avoid_cycling_back_to(const ProgressMeasures& current_progress) = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;

   protected:
      const double armijo_decrease_fraction; /*!< Sufficient reduction constant */
      const double armijo_tolerance;
      const bool protect_actual_reduction_against_roundoff;
      const double protected_actual_reduction_macheps_coefficient;
      const double small_infeasibility_factor;

      [[nodiscard]] bool armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const;
      [[nodiscard]] double compute_actual_reduction(double current_measure, double trial_measure) const;
   };
} // namespace

#endif // UNO_GLOBALIZATIONSTRATEGY_H