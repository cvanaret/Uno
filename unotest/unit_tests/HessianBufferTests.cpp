// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression tests: a Lagrangian Hessian with non-finite values (e.g. the second derivative of |x|^1.852 at 0, as in
// the CUTEst gas network NET3) must not reach the KKT matrix, where it corrupts the factorization and makes the inertia
// correction fail. It is handled like an evaluation failure: the previous Hessian (zero at the first evaluation) is kept.

#include <cmath>
#include <limits>
#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "ingredients/hessian_models/ExactHessian.hpp"

using namespace uno;

namespace {
   // min 0.5 * a * (x0^2 + x1^2) (lower triangle of the Hessian: (0,0), (1,0), (1,1)), where a is controlled by the test.
   // A NaN or infinite a makes the Hessian non-finite
   class ControllableHessianModel: public DenseTestModel {
   public:
      double& curvature;

      explicit ControllableHessianModel(double& curvature):
         DenseTestModel("controllable_hessian", 2, 0, ProblemType::NONLINEAR,
            [&curvature](const Vector<double>& x) { return 0.5 * curvature * (x[0]*x[0] + x[1]*x[1]); },
            [&curvature](const Vector<double>& x, Vector<double>& g) { g[0] = curvature * x[0]; g[1] = curvature * x[1]; },
            [](const Vector<double>&, Vector<double>&) { },
            [](const Vector<double>&, double*) { },
            [&curvature](const Vector<double>&, double rho, const Vector<double>&, double* H) {
               H[0] = rho * curvature; H[1] = 0.; H[2] = rho * curvature;
            },
            {-Inf, -Inf}, {Inf, Inf}, {}, {}, {1., 1.}),
         curvature(curvature) { }
   };

   void evaluate(HessianBuffer& buffer, const Model& model, Vector<double>& hessian_values) {
      const Vector<double> x{1., 1.};
      const Vector<double> multipliers(0);
      buffer.evaluate_hessian(model, x, 1., multipliers, hessian_values.view());
   }
} // namespace

TEST(HessianBuffer, FiniteHessianIsCopied) {
   double curvature = 3.;
   const ControllableHessianModel model(curvature);
   HessianBuffer buffer(model);
   Vector<double> hessian_values(3, -1.);
   evaluate(buffer, model, hessian_values);
   EXPECT_EQ(hessian_values[0], 3.);
   EXPECT_EQ(hessian_values[1], 0.);
   EXPECT_EQ(hessian_values[2], 3.);
}

TEST(HessianBuffer, NonFiniteInitialHessianIsReplacedByZero) {
   double curvature = std::numeric_limits<double>::quiet_NaN();
   const ControllableHessianModel model(curvature);
   HessianBuffer buffer(model);
   Vector<double> hessian_values(3, -1.);
   evaluate(buffer, model, hessian_values);
   for (size_t index = 0; index < 3; ++index) {
      EXPECT_EQ(hessian_values[index], 0.);
   }
}

TEST(HessianBuffer, NonFiniteHessianKeepsThePreviousOne) {
   double curvature = 2.;
   const ControllableHessianModel model(curvature);
   HessianBuffer buffer(model);
   Vector<double> hessian_values(3, 0.);
   evaluate(buffer, model, hessian_values);
   ASSERT_EQ(hessian_values[0], 2.);

   for (const double invalid_curvature: {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
      curvature = invalid_curvature;
      evaluate(buffer, model, hessian_values);
      EXPECT_EQ(hessian_values[0], 2.);
      EXPECT_EQ(hessian_values[1], 0.);
      EXPECT_EQ(hessian_values[2], 2.);
   }
}
