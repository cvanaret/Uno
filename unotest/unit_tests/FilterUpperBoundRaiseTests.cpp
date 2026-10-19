// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Adaptive infeasibility upper bound of the Waechter filter method: the upper bound is raised after a streak of line
// searches whose step length was limited by the upper bound

#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/globalization_strategies/switching_methods/filter_methods/WaechterFilterMethod.hpp"
#include "optimization/Iterate.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   Options make_options(size_t trigger, size_t max_raises = 4) {
      Options options;
      DefaultOptions::load(options);
      options.set_double("filter_ubd", 1e4);
      options.set_double("filter_fact", 1e4);
      options.set_integer("filter_ubd_raise_trigger", static_cast<uno_int>(trigger));
      options.set_double("filter_ubd_raise_factor", 100.);
      options.set_integer("filter_ubd_max_raises", static_cast<uno_int>(max_raises));
      return options;
   }

   ProgressMeasures progress(double infeasibility, double objective) {
      return {infeasibility, [=](double) { return objective; }, 0., 1.};
   }

   struct Fixture {
      Statistics statistics{};
      WaechterFilterMethod method;
      double objective{0.};

      explicit Fixture(const Options& options): method(options) {
         this->statistics.add_column("Status", Statistics::string_width, 0);
         Iterate initial_iterate(1, 1);
         initial_iterate.progress = progress(0., 0.); // feasible initial point: upper bound = 1e4
         this->method.initialize(this->statistics, initial_iterate);
      }

      bool test(const ProgressMeasures& current, double trial_infeasibility, double trial_objective) {
         // zero predicted reduction: the switching condition does not hold (h-type)
         return this->method.is_iterate_acceptable(this->statistics, current, progress(trial_infeasibility, trial_objective),
            progress(0., 0.), 1.);
      }

      // a line search whose full step is rejected by the upper bound and whose shorter step is accepted
      void upper_bound_limited_line_search() {
         const ProgressMeasures current = progress(1., this->objective);
         ASSERT_FALSE(this->test(current, 2e4, this->objective - 1.));
         ASSERT_TRUE(this->test(current, 0.5, this->objective - 1.));
         this->objective -= 1.;
      }
   };
} // namespace

TEST(FilterUpperBoundRaise, RaisedAfterStreak) {
   Fixture fixture(make_options(3));
   for (size_t line_search = 0; line_search < 3; ++line_search) {
      fixture.upper_bound_limited_line_search();
   }
   // the streak is closed when the next line search starts: its full step now passes the raised upper bound (1e6)
   EXPECT_TRUE(fixture.test(progress(1., fixture.objective), 2e4, fixture.objective - 1e3));
}

TEST(FilterUpperBoundRaise, NoRaiseBeforeStreak) {
   Fixture fixture(make_options(3));
   for (size_t line_search = 0; line_search < 2; ++line_search) {
      fixture.upper_bound_limited_line_search();
   }
   EXPECT_FALSE(fixture.test(progress(1., fixture.objective), 2e4, fixture.objective - 1e3));
}

TEST(FilterUpperBoundRaise, Disabled) {
   Fixture fixture(make_options(0));
   for (size_t line_search = 0; line_search < 10; ++line_search) {
      fixture.upper_bound_limited_line_search();
   }
   EXPECT_FALSE(fixture.test(progress(1., fixture.objective), 2e4, fixture.objective - 1e3));
}

// a line search with a rejection by another test (here: the current iterate) is not limited by the upper bound
TEST(FilterUpperBoundRaise, OtherRejectionBreaksStreak) {
   Fixture fixture(make_options(3));
   fixture.upper_bound_limited_line_search();
   fixture.upper_bound_limited_line_search();
   {
      const ProgressMeasures current = progress(1., fixture.objective);
      ASSERT_FALSE(fixture.test(current, 2e4, fixture.objective - 1.)); // upper bound
      ASSERT_FALSE(fixture.test(current, 1., fixture.objective + 1.)); // not acceptable wrt current iterate
      ASSERT_TRUE(fixture.test(current, 0.5, fixture.objective - 1.));
      fixture.objective -= 1.;
   }
   fixture.upper_bound_limited_line_search();
   EXPECT_FALSE(fixture.test(progress(1., fixture.objective), 2e4, fixture.objective - 1e3));
}

TEST(FilterUpperBoundRaise, MaximumNumberOfRaises) {
   Fixture fixture(make_options(1, 1));
   fixture.upper_bound_limited_line_search(); // closed by the next line search: 1e4 -> 1e6
   // full step 2e6 rejected even at 1e6, shorter step accepted: limited, but the maximum number of raises is reached
   for (size_t line_search = 0; line_search < 3; ++line_search) {
      const ProgressMeasures current = progress(1., fixture.objective);
      ASSERT_FALSE(fixture.test(current, 2e6, fixture.objective - 1.));
      ASSERT_TRUE(fixture.test(current, 0.5, fixture.objective - 1.));
      fixture.objective -= 1.;
   }
   EXPECT_FALSE(fixture.test(progress(1., fixture.objective), 2e6, fixture.objective - 1e5));
   EXPECT_TRUE(fixture.test(progress(1., fixture.objective), 2e4, fixture.objective - 1e3));
}
