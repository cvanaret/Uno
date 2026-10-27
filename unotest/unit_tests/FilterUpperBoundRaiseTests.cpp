// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Adaptive infeasibility upper bound of the Waechter filter method: the upper bound is raised after a streak of steps
// limited by the upper bound only

#include <gtest/gtest.h>
#include <initializer_list>
#include <utility>
#include "ingredients/globalization_mechanisms/RejectionHistory.hpp"
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

   struct Scenario {
      Statistics statistics{};
      WaechterFilterMethod method;
      double objective{0.};

      explicit Scenario(const Options& options): method(options) {
         this->statistics.add_column("Status", Statistics::string_width, 0);
         Iterate initial_iterate(1, 1);
         initial_iterate.progress = progress(0., 0.); // feasible initial point: upper bound = 1e4
         this->method.initialize(this->statistics, initial_iterate);
      }

      Decision test(double trial_infeasibility, double trial_objective) {
         // zero predicted reduction: the switching condition does not hold (h-type)
         return this->method.is_iterate_acceptable(this->statistics, progress(1., this->objective),
            progress(trial_infeasibility, trial_objective), progress(0., 0.), 1.);
      }

      // emulates a backtracking line search: the trial iterates (infeasibility, objective decrease) are tested in turn
      // until one is accepted, then the filter is notified. Returns whether a trial iterate was accepted
      bool step(std::initializer_list<std::pair<double, double>> trials) {
         RejectionHistory rejection_history;
         bool accepted = false;
         for (const auto& [trial_infeasibility, objective_decrease]: trials) {
            const Decision decision = this->test(trial_infeasibility, this->objective - objective_decrease);
            rejection_history.record(decision);
            if (decision.accepted) {
               this->objective -= objective_decrease;
               accepted = true;
               break;
            }
         }
         this->method.notify_step_limited_by(rejection_history.sole_reason());
         return accepted;
      }

      // the full step is rejected by the upper bound and a shorter step is accepted
      void upper_bound_limited_step(double full_step_infeasibility = 2e4) {
         ASSERT_TRUE(this->step({{full_step_infeasibility, 1.}, {0.5, 1.}}));
      }

      [[nodiscard]] bool full_step_acceptable(double infeasibility) {
         return this->test(infeasibility, this->objective - 1e3).accepted;
      }
   };
} // namespace

TEST(RejectionHistory, SoleReason) {
   RejectionHistory empty;
   EXPECT_EQ(empty.sole_reason(), RejectionReason::NONE);

   RejectionHistory accepted;
   accepted.record(Decision::accept());
   EXPECT_EQ(accepted.sole_reason(), RejectionReason::NONE);

   RejectionHistory single;
   single.record(Decision::reject(RejectionReason::INFEASIBILITY_UPPER_BOUND));
   single.record(Decision::reject(RejectionReason::INFEASIBILITY_UPPER_BOUND));
   single.record(Decision::accept());
   EXPECT_EQ(single.sole_reason(), RejectionReason::INFEASIBILITY_UPPER_BOUND);

   RejectionHistory mixed;
   mixed.record(Decision::reject(RejectionReason::INFEASIBILITY_UPPER_BOUND));
   mixed.record(Decision::reject(RejectionReason::FILTER));
   mixed.record(Decision::reject(RejectionReason::INFEASIBILITY_UPPER_BOUND));
   EXPECT_EQ(mixed.sole_reason(), RejectionReason::NONE);
}

TEST(FilterUpperBoundRaise, RaisedAfterStreak) {
   Scenario scenario(make_options(3));
   for (size_t step = 0; step < 3; ++step) {
      scenario.upper_bound_limited_step();
   }
   // the full step now passes the raised upper bound (1e6)
   EXPECT_TRUE(scenario.full_step_acceptable(2e4));
}

TEST(FilterUpperBoundRaise, NoRaiseBeforeStreak) {
   Scenario scenario(make_options(3));
   for (size_t step = 0; step < 2; ++step) {
      scenario.upper_bound_limited_step();
   }
   EXPECT_FALSE(scenario.full_step_acceptable(2e4));
}

TEST(FilterUpperBoundRaise, Disabled) {
   Scenario scenario(make_options(0));
   for (size_t step = 0; step < 10; ++step) {
      scenario.upper_bound_limited_step();
   }
   EXPECT_FALSE(scenario.full_step_acceptable(2e4));
}

// a step with a rejection by another test (here: the current iterate) is not limited by the upper bound only
TEST(FilterUpperBoundRaise, OtherRejectionBreaksStreak) {
   Scenario scenario(make_options(3));
   scenario.upper_bound_limited_step();
   scenario.upper_bound_limited_step();
   // upper bound, then not acceptable wrt the current iterate, then accepted
   ASSERT_TRUE(scenario.step({{2e4, 1.}, {1., -1.}, {0.5, 1.}}));
   scenario.upper_bound_limited_step();
   EXPECT_FALSE(scenario.full_step_acceptable(2e4));
}

// a full step accepted right away is not limited by the upper bound
TEST(FilterUpperBoundRaise, FullStepBreaksStreak) {
   Scenario scenario(make_options(3));
   scenario.upper_bound_limited_step();
   scenario.upper_bound_limited_step();
   ASSERT_TRUE(scenario.step({{0.5, 1.}}));
   scenario.upper_bound_limited_step();
   EXPECT_FALSE(scenario.full_step_acceptable(2e4));
}

// a failed step whose trial iterates were all rejected by the upper bound is limited by the upper bound only
TEST(FilterUpperBoundRaise, FailedStepCounts) {
   Scenario scenario(make_options(2));
   scenario.upper_bound_limited_step();
   ASSERT_FALSE(scenario.step({{2e4, 1.}, {1.5e4, 1.}}));
   EXPECT_TRUE(scenario.full_step_acceptable(2e4));
}

TEST(FilterUpperBoundRaise, MaximumNumberOfRaises) {
   Scenario scenario(make_options(1, 1));
   scenario.upper_bound_limited_step(); // 1e4 -> 1e6
   // full step 2e6 rejected even at 1e6, shorter step accepted: limited, but the maximum number of raises is reached
   for (size_t step = 0; step < 3; ++step) {
      scenario.upper_bound_limited_step(2e6);
   }
   EXPECT_FALSE(scenario.full_step_acceptable(2e6));
   EXPECT_TRUE(scenario.full_step_acceptable(2e4));
}
