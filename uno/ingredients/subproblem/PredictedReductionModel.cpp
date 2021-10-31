#include "PredictedReductionModel.hpp"

// PredictedReductionModel
PredictedReductionModel::PredictedReductionModel(double full_step_value,
      const std::function<std::function<double (double step_length)> ()>& partial_step_precomputation):
      full_step_predicted_reduction(full_step_value), partial_step_precomputation(partial_step_precomputation) {
}

double PredictedReductionModel::evaluate(double step_length) {
   if (step_length == 1.) {
      return this->full_step_predicted_reduction;
   }
   else {
      // unique evaluation of partial_step_precomputation
      if (this->partial_step_predicted_reduction == nullptr) {
         // precompute the expensive stuff
         this->partial_step_predicted_reduction = this->partial_step_precomputation();
      }
      return this->partial_step_predicted_reduction(step_length);
   }
}