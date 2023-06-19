// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NONMONOTONEFILTER_H
#define UNO_NONMONOTONEFILTER_H

#include "Filter.hpp"

class NonmonotoneFilter : public Filter {
public:
   explicit NonmonotoneFilter(const Options& options);

   void add(double infeasibility_measure, double optimality_measure) override;
   bool acceptable(double infeasibility_measure, double optimality_measure) override;
   bool acceptable_wrt_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
         double trial_optimality_measure) override;
   double compute_actual_reduction(double current_optimality_measure, double current_infeasibility_measure, double trial_optimality_measure) override;

protected:
   const size_t max_number_dominated_entries; /*!< Memory of filter */

   size_t compute_number_dominated_entries(double infeasibility_measure, double optimality_measure);
};

#endif // UNO_NONMONOTONEFILTER_H