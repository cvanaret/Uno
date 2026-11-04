// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// selection of the inverse penalty parameter τ of the merit function τ (objective + auxiliary) + infeasibility:
// τ is the largest value (nonincreasing until reset) such that pred(τ) ≥ σ Δinfeasibility, with σ = 0.1 (default)

#include <memory>
#include <tuple>
#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/MeritFunction.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   // objective measure: objective_multiplier * objective
   ProgressMeasures progress(double infeasibility, double objective, double auxiliary = 0.) {
      return {infeasibility, [=](double objective_multiplier) { return objective_multiplier * objective; }, auxiliary, 1.};
   }

   // exposes the inverse penalty parameter
   class TestMeritFunction: public MeritFunction {
   public:
      using MeritFunction::MeritFunction;
      [[nodiscard]] double get_inverse_penalty_parameter() const { return this->inverse_penalty_parameter; }
   };

   class MeritFunctionTest: public ::testing::Test {
   protected:
      Options options;
      Statistics statistics{};
      std::unique_ptr<TestMeritFunction> merit_function;

      void SetUp() override {
         DefaultOptions::load(this->options);
         this->statistics.add_column("Status", Statistics::string_width, 0);
         this->merit_function = std::make_unique<TestMeritFunction>(this->options);
         Iterate initial_iterate(1, 1);
         initial_iterate.progress = progress(1., 0.);
         this->merit_function->initialize(this->statistics, initial_iterate);
      }

      Decision test(const ProgressMeasures& current, const ProgressMeasures& trial, const ProgressMeasures& predicted,
            double objective_multiplier = 1.) {
         return this->merit_function->is_iterate_acceptable(this->statistics, current, trial, predicted, objective_multiplier);
      }
   };
} // namespace

TEST_F(MeritFunctionTest, InitialValueIsOne) {
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 1.);
}

TEST_F(MeritFunctionTest, UnchangedWhenOptimalityMeasureDecreases) {
   EXPECT_TRUE(this->test(progress(1., 0.), progress(0., -0.5), progress(1., 0.5)).accepted);
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 1.);
}

TEST_F(MeritFunctionTest, UnchangedWhenSufficientDecreaseHolds) {
   // pred(1) = -0.5 + 1 = 0.5 ≥ 0.1
   std::ignore = this->test(progress(1., 0.), progress(0., 0.5), progress(1., -0.5));
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 1.);
}

TEST_F(MeritFunctionTest, DecreasedToSatisfySufficientDecrease) {
   // τ = (1 - 0.1) * 1 / 9 = 0.1. With τ = 1, the merit function would increase (1 -> 5) and the step be rejected;
   // with τ = 0.1: merit 1 -> 0.5, pred(τ) = 0.1 (-9) + 1 = 0.1
   EXPECT_TRUE(this->test(progress(1., 0.), progress(0., 5.), progress(1., -9.)).accepted);
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 0.1);
}

TEST_F(MeritFunctionTest, UnchangedWhenInfeasibilityDoesNotDecrease) {
   std::ignore = this->test(progress(0., 0.), progress(0., 1.), progress(0., -9.));
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 1.);
}

TEST_F(MeritFunctionTest, AuxiliaryMeasureIsScaled) {
   // the auxiliary reduction (e.g. barrier terms) enters the optimality measure: τ = 0.9 / 9 = 0.1
   EXPECT_TRUE(this->test(progress(1., 0., 0.), progress(0., 0., 5.), progress(1., 0., -9.)).accepted);
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 0.1);
}

TEST_F(MeritFunctionTest, ObjectiveMeasureUsesObjectiveMultiplier) {
   // objective multiplier 0 (feasibility problem): the objective does not contribute, τ is unchanged
   std::ignore = this->test(progress(1., 0.), progress(0., 5.), progress(1., -9.), /* objective_multiplier = */ 0.);
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 1.);
}

TEST_F(MeritFunctionTest, IsMonotoneUntilReset) {
   std::ignore = this->test(progress(1., 0.), progress(0., 5.), progress(1., -9.));
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 0.1);
   // threshold 0.9 > 0.1: τ is not increased
   std::ignore = this->test(progress(1., 0.), progress(0., 0.5), progress(1., -1.));
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 0.1);
   this->merit_function->reset();
   EXPECT_DOUBLE_EQ(this->merit_function->get_inverse_penalty_parameter(), 1.);
}

TEST_F(MeritFunctionTest, InsufficientActualReductionIsRejected) {
   // τ = 0.1: merit 1 -> 0.1 * 50 + 0 = 5
   EXPECT_FALSE(this->test(progress(1., 0.), progress(0., 50.), progress(1., -9.)).accepted);
}
