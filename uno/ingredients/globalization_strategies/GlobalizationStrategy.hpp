// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONSTRATEGY_H
#define UNO_GLOBALIZATIONSTRATEGY_H

#include <string>
#include "tools/Infinity.hpp"
#include "Decision.hpp"

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
      [[nodiscard]] virtual Decision is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reductions, double objective_multiplier) = 0;
      // called by the globalization mechanism once a step is complete: sole_rejection_reason is the reason shared by all
      // the rejected trial iterates of the step (NONE if no trial iterate was rejected or if the reasons differ)
      virtual void notify_step_limited_by(RejectionReason /*sole_rejection_reason*/) { }
      [[nodiscard]] virtual bool is_infeasibility_sufficiently_reduced(const Iterate& trial_iterate,
         double reference_infeasibility) const = 0;

      virtual void reset() = 0;

      // theta_min: below this infeasibility, the switching rule may declare f-type iterations. The line search uses the
      // same threshold to compute its minimum step length (IPOPT, Wächter & Biegler 2006, eq. (23))
      [[nodiscard]] virtual double small_infeasibility_threshold() const;

      virtual void avoid_cycling_back_to(const ProgressMeasures& current_progress) = 0;

      // set by the globalization mechanism before testing a trial iterate (trust region: h-type iterates are allowed
      // only if the step is short wrt the current infeasibility). Strategies without h-type iterates ignore it
      void set_h_type_iterates_allowed(bool allowed) { this->h_type_iterates_allowed = allowed; }
      // reference infeasibility tau of FLT's Lemma 4 (smallest infeasibility of the filter, width of the funnel): the
      // trust region bounds the length of h-type steps by factor tau^nu. Infinite for strategies without h-type iterates
      [[nodiscard]] virtual double h_type_reference_infeasibility() const { return Inf; }

      [[nodiscard]] virtual std::string get_name() const = 0;

   protected:
      const double armijo_decrease_fraction; /*!< Sufficient reduction constant */
      const double armijo_tolerance;
      const bool protect_actual_reduction_against_roundoff;
      const double protected_actual_reduction_macheps_coefficient;
      const double small_infeasibility_factor;
      bool h_type_iterates_allowed{true};

      [[nodiscard]] bool armijo_sufficient_decrease(double predicted_reduction, double actual_reduction) const;
      [[nodiscard]] double compute_actual_reduction(double current_measure, double trial_measure) const;
   };
} // namespace

#endif // UNO_GLOBALIZATIONSTRATEGY_H