// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FLETCHERFILTERMETHOD_H
#define UNO_FLETCHERFILTERMETHOD_H

#include "FilterMethod.hpp"

namespace uno {
   class FletcherFilterMethod : public FilterMethod {
   public:
      explicit FletcherFilterMethod(const Options& options);
      ~FletcherFilterMethod() override;

      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) override;
      [[nodiscard]] bool is_infeasibility_sufficiently_reduced(const ProgressMeasures& reference_progress,
         const ProgressMeasures& trial_progress) const override;

      [[nodiscard]] std::string get_name() const override;
   };
} // namespace

#endif // UNO_FLETCHERFILTERMETHOD_H