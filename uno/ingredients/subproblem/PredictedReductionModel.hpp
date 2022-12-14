// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREDICTEDREDUCTIONMODEL_H
#define UNO_PREDICTEDREDUCTIONMODEL_H

#include <functional>

class PredictedReductionModel {
public:
   explicit PredictedReductionModel(std::function<double(double step_length)> step_precomputation);
   double evaluate(double step_length);

private:
   // predicted reduction, function of the step length
   std::function<double (double step_length)> predicted_reduction;
};

#endif // UNO_PREDICTEDREDUCTIONMODEL_H
