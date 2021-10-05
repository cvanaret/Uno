#ifndef PREDICTEDREDUCTIONMODEL_H
#define PREDICTEDREDUCTIONMODEL_H

#include <functional>

struct PredictedReductionModel {
   // predicted reduction for a full step
   const double full_step_predicted_reduction;

   // this function, when evaluated, precomputes expensive quantities and returns a function of the step length
   const std::function<std::function<double (double step_length)> ()> partial_step_precomputation;

   // predicted reduction, function of the step length
   std::function<double (double step_length)> partial_step_predicted_reduction{nullptr};

   PredictedReductionModel(double full_step_value, const std::function<std::function<double(double step_length)>()>& partial_step_precomputation);
   double evaluate(double step_length);
};

#endif // PREDICTEDREDUCTIONMODEL_H
