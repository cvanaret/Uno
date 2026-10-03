// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <utility>
#include <vector>
#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/switching_methods/filter_methods/Filter.hpp"
#include "options/Options.hpp"

using namespace uno;

namespace {
   Options make_filter_options(double beta, double gamma, uno_int capacity) {
      Options options;
      options.set_double("filter_beta", beta);
      options.set_double("filter_gamma", gamma);
      options.set_integer("filter_capacity", capacity);
      options.set_bool("filter_margin_uses_entry_infeasibility", true); // a la IPOPT
      return options;
   }

   // exposes the stored (infeasibility, objective) entries, in storage order
   class InspectableFilter: public Filter {
   public:
      using Filter::Filter;

      [[nodiscard]] std::vector<std::pair<double, double>> entries() const {
         std::vector<std::pair<double, double>> result;
         for (size_t index = 0; index < this->number_entries; ++index) {
            result.emplace_back(this->infeasibility[index], this->objective[index]);
         }
         return result;
      }
   };

   // beta close to 1 (IPOPT preset) makes the margin tiny, which is where the insertion bug shows up
   constexpr double filter_beta = 0.99999;
   constexpr double filter_gamma = 1e-5;
   constexpr uno_int filter_capacity = 50;

   // the new entry has a smaller infeasibility than the existing one (0.999995 < 1.0), but lies within
   // the beta margin (0.999995 > beta * 1.0 = 0.99999). It does not dominate (1.0, 0) nor is dominated by it
   constexpr double current_infeasibility = 1.0;
   constexpr double current_objective = 0.;
   constexpr double trial_infeasibility = 0.999995;
   constexpr double trial_objective = 5.;
}

// the stored entries must be sorted by increasing infeasibility and decreasing objective (raw values, no margin)
TEST(Filter, EntryWithinMarginIsInsertedInSortedOrder) {
   InspectableFilter filter(make_filter_options(filter_beta, filter_gamma, filter_capacity));
   filter.add(current_infeasibility, current_objective);
   filter.add(trial_infeasibility, trial_objective);

   const std::vector<std::pair<double, double>> expected_entries{
      {trial_infeasibility, trial_objective},
      {current_infeasibility, current_objective}
   };
   EXPECT_EQ(filter.entries(), expected_entries);
}

// the left-most entry must be the one with the smallest infeasibility
TEST(Filter, SmallestInfeasibilityAfterInsertionWithinMargin) {
   Filter filter(make_filter_options(filter_beta, filter_gamma, filter_capacity));
   filter.add(current_infeasibility, current_objective);
   filter.add(trial_infeasibility, trial_objective);

   EXPECT_DOUBLE_EQ(filter.get_smallest_infeasibility(), trial_infeasibility);
}

// (1.5, 3) is dominated by the entry (1.0, 0) and must be rejected. With the entries out of order, the
// acceptability test compares the trial objective against the wrong entry (objective 5) and accepts it
TEST(Filter, DominatedTrialIsRejectedAfterInsertionWithinMargin) {
   Filter filter(make_filter_options(filter_beta, filter_gamma, filter_capacity));
   filter.add(current_infeasibility, current_objective);
   filter.add(trial_infeasibility, trial_objective);

   EXPECT_FALSE(filter.filter_acceptable(1.5, 3.));
}