// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WAECHTERFILTERMETHOD_H
#define UNO_WAECHTERFILTERMETHOD_H

#include "FilterMethod.hpp"

class WaechterFilterMethod : public FilterMethod {
public:
   explicit WaechterFilterMethod(Statistics& statistics, const Options& options);

   void initialize(const Iterate& initial_iterate) override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const Iterate& trial_iterate, const ProgressMeasures& current_progress_measures,
         const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction, double objective_multiplier) override;

protected:
   double initial_infeasibility{INF<double>};
};

#endif // UNO_WAECHTERFILTERMETHOD_H