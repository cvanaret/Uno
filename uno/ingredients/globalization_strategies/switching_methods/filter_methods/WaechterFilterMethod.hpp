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
      ~WaechterFilterMethod();

      void initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) override;
      [[nodiscard]] bool is_infeasibility_sufficiently_reduced(const ProgressMeasures& reference_progress, const ProgressMeasures& trial_progress) const override;

   protected:
      double initial_infeasibility{INF<double>};
      const double sufficient_infeasibility_decrease_factor;

      [[nodiscard]] bool is_regular_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
            const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) override;
   };
} // namespace

#endif // UNO_WAECHTERFILTERMETHOD_H
