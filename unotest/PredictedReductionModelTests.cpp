#include <gtest/gtest.h>
#include "ingredients/subproblem/Subproblem.hpp"

TEST(PredictedReductionModel, Correctness) {
   const double full_step_predicted_reduction = 1.;
   PredictedReductionModel predicted_reduction_model(full_step_predicted_reduction, [&]() {
      // expensive computations here
      const double quadratic_term = 2.;
      const double linear_term = -3.;
      return [=](double step_length) {
         return -step_length*(linear_term + step_length*quadratic_term);
      };
   });
   ASSERT_EQ(predicted_reduction_model.evaluate(1.), 1.);
   ASSERT_EQ(predicted_reduction_model.evaluate(0.25), 0.625);
}

TEST(PredictedReductionModel, Laziness) {
   const double full_step_predicted_reduction = 1.;
   size_t counter_expensive_computations = 0;
   PredictedReductionModel predicted_reduction_model(full_step_predicted_reduction, [&]() {
      // expensive computations here
      const double quadratic_term = 2.;
      const double linear_term = -3.;
      counter_expensive_computations++;
      return [=](double step_length) {
         return -step_length*(linear_term + step_length*quadratic_term);
      };
   });
   predicted_reduction_model.evaluate(1.);
   ASSERT_EQ(counter_expensive_computations, 0);
}

TEST(PredictedReductionModel, Precomputation) {
   const double full_step_predicted_reduction = 1.;
   size_t counter_expensive_computations = 0;
   PredictedReductionModel predicted_reduction_model(full_step_predicted_reduction, [&]() {
      // expensive computations here
      const double quadratic_term = 2.;
      const double linear_term = -3.;
      counter_expensive_computations++;
      return [=](double step_length) {
         return -step_length*(linear_term + step_length*quadratic_term);
      };
   });
   predicted_reduction_model.evaluate(1.);
   predicted_reduction_model.evaluate(0.5);
   predicted_reduction_model.evaluate(0.25);
   ASSERT_EQ(counter_expensive_computations, 1);
}