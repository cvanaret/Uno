// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREDICTEDREDUCTIONMODEL_H
#define UNO_PREDICTEDREDUCTIONMODEL_H

#include <functional>

struct PredictedReduction {
   double infeasibility;
   std::function<double(double objective_multiplier)> scaled_optimality;
   double unscaled_optimality;
};

#endif // UNO_PREDICTEDREDUCTIONMODEL_H
