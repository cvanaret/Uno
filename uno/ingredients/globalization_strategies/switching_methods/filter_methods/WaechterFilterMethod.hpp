// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WAECHTERFILTERMETHOD_H
#define UNO_WAECHTERFILTERMETHOD_H

#include "FilterMethod.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   class WaechterFilterMethod : public FilterMethod {
   public:
      explicit WaechterFilterMethod(const Options& options);
      ~WaechterFilterMethod() override = default;

      void initialize(Statistics& statistics, const Iterate& initial_iterate) override;
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reductions, double objective_multiplier) override;
      [[nodiscard]] bool is_infeasibility_sufficiently_reduced(const Iterate& trial_iterate,
         double reference_infeasibility) const override;
      [[nodiscard]] double small_infeasibility_threshold() const override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      double initial_infeasibility{Inf};
      const double sufficient_infeasibility_decrease_factor;
      bool last_rejection_due_to_filter{false};
      // filter resets
      size_t number_successive_filter_rejections{0};
      const size_t filter_reset_iteration_threshold; // >= 1
      size_t number_filter_resets{0};
      const size_t max_number_filter_resets; // >= 0

      // adaptive infeasibility upper bound (POUNCE, pounce#546): after a streak of line searches whose step length was
      // limited by the upper bound, the upper bound is raised
      const size_t upper_bound_raise_trigger; // 0: disabled
      const double upper_bound_raise_factor; // > 1
      const size_t max_number_upper_bound_raises;
      size_t number_upper_bound_raises{0};
      size_t upper_bound_limited_line_search_streak{0};
      // bookkeeping of the current line search. A new line search is detected when the current progress measures change
      bool line_search_started{false};
      double line_search_current_infeasibility{0.};
      double line_search_current_merit{0.};
      size_t number_trials_this_line_search{0};
      size_t number_rejected_trials_this_line_search{0};
      size_t number_upper_bound_rejections_this_line_search{0};
      bool first_trial_rejected_by_upper_bound{false};

      void register_trial(double current_infeasibility, double current_merit);
      void end_of_line_search();
   };
} // namespace

#endif // UNO_WAECHTERFILTERMETHOD_H
