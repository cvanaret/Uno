// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Exit test of the restoration phase in the Waechter filter method: the restoration iterate is handed back to the
// optimality phase if
// - its primal infeasibility is sufficiently reduced: inf_pr <= max(kappa_resto inf_pr_ref, min(tol, constr_viol_tol))
//   (IPOPT, RestoConvergenceCheck);
// - it is acceptable to the optimality filter (IPOPT, RestoFilterConvergenceCheck);
// - its infeasibility does not exceed theta_max: in the theory (Wächter & Biegler), the initial filter is
//   F_0 = {theta >= theta_max} and the new iterate must lie outside the filter. IPOPT does not test theta_max there

#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/globalization_strategies/switching_methods/filter_methods/WaechterFilterMethod.hpp"
#include "optimization/Iterate.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   ProgressMeasures progress(double infeasibility, double objective) {
      return {infeasibility, [=](double) { return objective; }, 0., 1.};
   }

   Iterate iterate(double infeasibility, double primal_infeasibility, double objective) {
      Iterate iterate(1, 1);
      iterate.progress = progress(infeasibility, objective);
      iterate.scaled_residuals.primal_infeasibility = primal_infeasibility;
      return iterate;
   }

   Options ipopt_options() {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      return options;
   }

   // optimality filter initialized at a point with infeasibility 1
   struct OptimalityFilter {
      OptimalityFilter() {
         this->statistics.add_column("Status", Statistics::string_width, 0);
         this->method.initialize(this->statistics, iterate(1., 1., 0.));
      }
      const Options options{ipopt_options()};
      Statistics statistics{};
      WaechterFilterMethod method{options};
   };
} // namespace

// the restoration phase was entered at an almost feasible point (inf_pr_ref = 1e-10, below primal_tolerance = 1e-8): a
// point that is feasible wrt the tolerance is handed back, although inf_pr > kappa_resto inf_pr_ref
TEST(WaechterRestorationExit, FeasiblePointWrtToleranceLeavesRestoration) {
   OptimalityFilter filter;
   const double reference_infeasibility = 1e-10;
   EXPECT_TRUE(filter.method.is_infeasibility_sufficiently_reduced(iterate(5e-9, 5e-9, 0.), reference_infeasibility));
   // above the tolerance: the sufficient reduction applies
   EXPECT_FALSE(filter.method.is_infeasibility_sufficiently_reduced(iterate(5e-8, 5e-8, 0.), reference_infeasibility));
}

// the primal infeasibility (max-norm) is sufficiently reduced, but the filter infeasibility (theta, in the progress norm)
// exceeds the upper bound theta_max of the optimality filter: the point lies in F_0 and stays in restoration
TEST(WaechterRestorationExit, PointAboveInfeasibilityUpperBoundStaysInRestoration) {
   OptimalityFilter filter;
   const double reference_infeasibility = 1.;
   EXPECT_TRUE(filter.method.is_infeasibility_sufficiently_reduced(iterate(0.5, 0.5, 0.), reference_infeasibility));
   EXPECT_FALSE(filter.method.is_infeasibility_sufficiently_reduced(iterate(1e20, 0.5, 0.), reference_infeasibility));
}
