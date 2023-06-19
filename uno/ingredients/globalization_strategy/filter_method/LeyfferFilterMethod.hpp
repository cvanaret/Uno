// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LEYFFERFILTERMETHOD_H
#define UNO_LEYFFERFILTERMETHOD_H

#include "FilterMethod.hpp"

class LeyfferFilterMethod : public FilterMethod {
public:
   LeyfferFilterMethod(Statistics& statistics, bool accept_when_switching_violated, const Options& options);

   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const Iterate& trial_iterate, const ProgressMeasures& current_progress_measures,
         const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction, double objective_multiplier) override;

protected:
   const bool accept_when_switching_violated;
};

#endif // UNO_LEYFFERFILTERMETHOD_H