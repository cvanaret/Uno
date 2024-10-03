// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LEYFFERFILTERMETHOD_H
#define UNO_LEYFFERFILTERMETHOD_H

#include "FilterMethod.hpp"

namespace uno {
   class FletcherFilterMethod : public FilterMethod {
   public:
      explicit FletcherFilterMethod(const Options& options);
      ~FletcherFilterMethod();

      [[nodiscard]] bool is_infeasibility_sufficiently_reduced(const ProgressMeasures& reference_progress, const ProgressMeasures& trial_progress) const override;

   protected:
      [[nodiscard]] bool is_regular_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
            const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) override;
   };
} // namespace

#endif // UNO_LEYFFERFILTERMETHOD_H
