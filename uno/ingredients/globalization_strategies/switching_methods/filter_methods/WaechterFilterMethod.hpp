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
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) override;
      [[nodiscard]] bool is_infeasibility_sufficiently_reduced(const ProgressMeasures& reference_progress,
         const ProgressMeasures& trial_progress) const override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      double initial_infeasibility{INF<double>};
      const double sufficient_infeasibility_decrease_factor;
      const double small_infeasibility_factor;
      bool last_rejection_due_to_filter{false};
      // filter resets
      size_t number_successive_filter_rejections{0};
      const size_t filter_reset_iteration_threshold; // >= 1
      size_t number_filter_resets{0};
      const size_t max_number_filter_resets; // >= 0
   };
} // namespace

#endif // UNO_WAECHTERFILTERMETHOD_H
