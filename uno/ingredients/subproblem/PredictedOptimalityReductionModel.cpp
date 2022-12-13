// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PredictedOptimalityReductionModel.hpp"

#include <utility>

// PredictedReductionModel
PredictedOptimalityReductionModel::PredictedOptimalityReductionModel(
      std::function<double (double)>  partial_step_precomputation):
      predicted_reduction(std::move(partial_step_precomputation)) {
}

double PredictedOptimalityReductionModel::evaluate(double step_length) {
      return this->predicted_reduction(step_length);
}