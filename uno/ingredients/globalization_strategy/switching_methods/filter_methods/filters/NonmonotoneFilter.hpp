// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NONMONOTONEFILTER_H
#define UNO_NONMONOTONEFILTER_H

#include "Filter.hpp"

namespace uno {
   class NonmonotoneFilter : public Filter {
   public:
      explicit NonmonotoneFilter(const Options& options);

      void add(double current_infeasibility, double current_objective) override;
      [[nodiscard]] bool acceptable(double trial_infeasibility, double trial_objective) override;
      [[nodiscard]] bool acceptable_wrt_current_iterate(double current_infeasibility, double current_objective, double trial_infeasibility,
            double trial_objective) const override;
      [[nodiscard]] double compute_actual_objective_reduction(double current_objective, double current_infeasibility, double trial_objective) override;

   protected:
      const size_t max_number_dominated_entries; /*!< Memory of filter */

      [[nodiscard]] size_t compute_number_dominated_entries(double trial_infeasibility, double trial_objective) const;
   };
} // namespace

#endif // UNO_NONMONOTONEFILTER_H
