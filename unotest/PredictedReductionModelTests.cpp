// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "ingredients/subproblem/PredictedOptimalityReductionModel.hpp"

TEST(PredictedOptimalityReductionModel, Correctness) {
   // expensive computations here
   const double quadratic_term = 2.;
   const double linear_term = -3.;
   PredictedOptimalityReductionModel predicted_reduction_model([=](double step_length) {
      return -step_length*(linear_term + step_length*quadratic_term);
   });
   ASSERT_EQ(predicted_reduction_model.evaluate(1.), 1.);
   ASSERT_EQ(predicted_reduction_model.evaluate(0.25), 0.625);
}
