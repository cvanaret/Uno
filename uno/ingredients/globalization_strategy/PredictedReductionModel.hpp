// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREDICTEDREDUCTIONMODEL_H
#define UNO_PREDICTEDREDUCTIONMODEL_H

#include <functional>
#include <string>

struct PredictedReductionModel {
   std::function<double(double step_length)> function;
   std::string text;

   double operator()(double step_length) const {
      return this->function(step_length);
   }
};

struct PredictedReductionModels {
   PredictedReductionModel infeasibility;
   PredictedReductionModel scaled_optimality;
   PredictedReductionModel unscaled_optimality;
};

#endif // UNO_PREDICTEDREDUCTIONMODEL_H
