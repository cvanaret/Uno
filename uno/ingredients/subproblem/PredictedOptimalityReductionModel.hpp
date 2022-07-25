// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREDICTEDOPTIMALITYREDUCTIONMODEL_H
#define UNO_PREDICTEDOPTIMALITYREDUCTIONMODEL_H

#include <functional>

class PredictedOptimalityReductionModel {
public:
   PredictedOptimalityReductionModel(double full_step_value, std::function<std::function<double(double step_length)>()> partial_step_precomputation);
   double evaluate(double step_length);

private:
   // predicted reduction for a full step
   const double full_step_predicted_reduction;

   // predicted reduction, function of the step length
   std::function<double (double step_length)> partial_step_predicted_reduction{nullptr};

   // this function, when evaluated, precomputes expensive quantities and returns a function of the step length
   const std::function<std::function<double (double step_length)> ()> partial_step_precomputation;
};

#endif // UNO_PREDICTEDOPTIMALITYREDUCTIONMODEL_H
