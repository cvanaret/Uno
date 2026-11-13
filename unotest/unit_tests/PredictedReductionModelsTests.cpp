// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/PredictedReductionModels.hpp"

using namespace uno;

namespace {
   PredictedReductionModels make_models() {
      return {
         [](double step_length) { return 2. * step_length; },
         // first-order term 3α, curvature term α²
         [](double step_length, double objective_multiplier, bool use_curvature_information) {
            return objective_multiplier * 3. * step_length - (use_curvature_information ? step_length * step_length : 0.);
         },
         [](double step_length) { return 5. * step_length; }
      };
   }
} // namespace

// the predicted reductions remain valid after the models are destroyed
TEST(PredictedReductionModels, ReductionsOutliveModels) {
   const ProgressMeasures reductions = make_models()(0.5, /* use_curvature_information = */ true);
   EXPECT_DOUBLE_EQ(reductions.infeasibility, 1.);
   EXPECT_DOUBLE_EQ(reductions.objective(1.), 1.5 - 0.25);
   EXPECT_DOUBLE_EQ(reductions.auxiliary, 2.5);
   EXPECT_DOUBLE_EQ(reductions.step_length, 0.5);
}

TEST(PredictedReductionModels, FirstOrderReductions) {
   const ProgressMeasures reductions = make_models()(0.5, /* use_curvature_information = */ false);
   EXPECT_DOUBLE_EQ(reductions.objective(2.), 3.);
   EXPECT_DOUBLE_EQ(make_models().optimality_first_order_reduction(1.), 3. + 5.);
}
