// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "ingredients/globalization_mechanisms/TrustRegionStrategy.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"

using namespace uno;

// Bug: TR_increase_factor and TR_decrease_factor were validated, but not TR_aggressive_decrease_factor. With a factor
// <= 1, the radius never shrinks (or grows) after an unbounded subproblem, and the trust-region loop never terminates.
TEST(TrustRegionStrategy, AggressiveDecreaseFactorIsValidated) {
   // min x^2  s.t.  x >= 1
   const DenseTestModel model("tr", 1, 1, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return x[0]*x[0]; },
      [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*x[0]; },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]; },
      [](const Vector<double>&, double* J) { J[0] = 1.; },
      [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; },
      {-Inf}, {Inf}, {1.}, {Inf}, {0.});
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "filtersqp");
   std::vector<OptionOverride> option_overrides;
   EXPECT_NO_THROW(TrustRegionStrategy(model, options, option_overrides));
   options.set_double("TR_aggressive_decrease_factor", 1.);
   EXPECT_THROW(TrustRegionStrategy(model, options, option_overrides), std::runtime_error);
}