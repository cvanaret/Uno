// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PredictedReductionModel.hpp"

#include <utility>

// PredictedReductionModel
PredictedReductionModel::PredictedReductionModel(
      std::function<double (double)>  partial_step_precomputation):
      predicted_reduction(std::move(partial_step_precomputation)) {
}

double PredictedReductionModel::evaluate(double step_length) {
      return this->predicted_reduction(step_length);
}