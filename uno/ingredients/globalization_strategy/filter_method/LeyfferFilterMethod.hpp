// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LEYFFERFILTERMETHOD_H
#define UNO_LEYFFERFILTERMETHOD_H

#include "FilterMethod.hpp"

class LeyfferFilterMethod : public FilterMethod {
public:
   LeyfferFilterMethod(bool solving_feasibility_problem, const Options& options);

   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress_measures,
         const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction, double objective_multiplier) override;
   [[nodiscard]] bool is_infeasibility_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress) const override;

protected:
   const bool solving_feasibility_problem;
};

#endif // UNO_LEYFFERFILTERMETHOD_H
