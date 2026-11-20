// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// The function scaling is computed at the point where the optimization starts: the initial point projected onto the
// bounds (SQP methods) or pushed strictly inside the bounds (interior-point methods), not at the initial point of the user

#include <optional>
#include <gtest/gtest.h>
#include "BoundActiveModel.hpp"
#include "Uno.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   // bounds: x0 >= 0, -1 <= x1 <= 1. The initial point (-1, 3) violates both. Records the point at which the objective
   // gradient is first evaluated (by the function scaling)
   class FirstGradientRecordingModel: public test::BoundActiveModel {
   public:
      void initial_primal_point(Vector<double>& x) const override {
         x[0] = -1.;
         x[1] = 3.;
      }
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         if (!this->first_gradient_point.has_value()) {
            this->first_gradient_point = x;
         }
         test::BoundActiveModel::evaluate_objective_gradient(x, gradient);
      }

      mutable std::optional<Vector<double>> first_gradient_point{};
   };

   // no iteration: only the reformulations of the model (function scaling included) are performed
   Options create_options(const std::string& inequality_handling_method) {
      Options options;
      DefaultOptions::load(options);
      options.set_string("inequality_handling_method", inequality_handling_method);
      options.set_bool("use_function_scaling", true);
      options.set_integer("max_iterations", 0);
      return options;
   }
} // namespace

TEST(FunctionScalingPoint, SQPScalesAtProjectedInitialPoint) {
   const FirstGradientRecordingModel model;
   const Options options = create_options("inequality_constrained");
   Logger::level = SILENT;

   Uno uno;
   (void) uno.solve(model, options);
   ASSERT_TRUE(model.first_gradient_point.has_value());
   EXPECT_DOUBLE_EQ((*model.first_gradient_point)[0], 0.);
   EXPECT_DOUBLE_EQ((*model.first_gradient_point)[1], 1.);
}

TEST(FunctionScalingPoint, InteriorPointScalesAtPushedInitialPoint) {
   const FirstGradientRecordingModel model;
   Options options = create_options("interior_point");
   options.set_double("barrier_push_variable_to_interior_k1", 1e-2);
   options.set_double("barrier_push_variable_to_interior_k2", 1e-2);
   Logger::level = SILENT;

   Uno uno;
   (void) uno.solve(model, options);
   ASSERT_TRUE(model.first_gradient_point.has_value());
   // pushed by 1e-2 inside the (slightly relaxed) bounds
   EXPECT_NEAR((*model.first_gradient_point)[0], 1e-2, 1e-6);
   EXPECT_NEAR((*model.first_gradient_point)[1], 1. - 1e-2, 1e-6);
}
