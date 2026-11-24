// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <random>
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
         for (size_t index = 0; index < this->number_entries(); ++index) {
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

// Wächter & Biegler's filter (entry margin, IPOPT preset): the taboo region of an entry (h, f) is
// {(h', f'): h' >= beta*h and f' > f - gamma*h}. The filter must forbid the union of the taboo regions of all the pairs
// added so far (F_k is included in F_{k+1}), whatever entries it prunes. Comparing raw (h, f) pairs to prune entries
// is not enough: a pair with a larger infeasibility has a lower corner f - gamma*h
namespace {
   constexpr double entry_margin_gamma = 1e-3;

   bool in_entry_margin_taboo_region(const std::vector<std::pair<double, double>>& added_pairs, double beta, double gamma,
         double trial_infeasibility, double trial_objective) {
      for (const auto& [entry_infeasibility, entry_objective]: added_pairs) {
         if (beta * entry_infeasibility <= trial_infeasibility && entry_objective - gamma * entry_infeasibility < trial_objective) {
            return true;
         }
      }
      return false;
   }
} // namespace

// (10, 0) is dominated by (1, 0) on the raw values, but its taboo region reaches lower objectives: f' > -10 gamma
TEST(Filter, RawDominatedPairKeepsItsEntryMarginTabooRegion) {
   Filter filter(make_filter_options(filter_beta, entry_margin_gamma, 0));
   filter.set_infeasibility_upper_bound(1e30);
   filter.add(1., 0.);
   filter.add(10., 0.);
   const double trial_objective = -5. * entry_margin_gamma;
   ASSERT_TRUE(in_entry_margin_taboo_region({{1., 0.}, {10., 0.}}, filter_beta, entry_margin_gamma, 10., trial_objective));
   EXPECT_FALSE(filter.filter_acceptable(10., trial_objective));
}

// same pairs in the reverse order: (10, 0) is removed when (1, 0) is added
TEST(Filter, PrunedEntryKeepsItsEntryMarginTabooRegion) {
   Filter filter(make_filter_options(filter_beta, entry_margin_gamma, 0));
   filter.set_infeasibility_upper_bound(1e30);
   filter.add(10., 0.);
   filter.add(1., 0.);
   EXPECT_FALSE(filter.filter_acceptable(10., -5. * entry_margin_gamma));
}

// random sequences of pairs: no point of the union of the taboo regions is acceptable, and every point outside of it is
TEST(Filter, EntryMarginFilterIsTheUnionOfTabooRegions) {
   std::mt19937 generator(7);
   std::uniform_real_distribution<double> uniform(0., 1.);
   size_t number_wrongly_accepted = 0;
   size_t number_wrongly_rejected = 0;
   for (size_t run = 0; run < 200; ++run) {
      Filter filter(make_filter_options(filter_beta, entry_margin_gamma, 0));
      filter.set_infeasibility_upper_bound(1e30);
      std::vector<std::pair<double, double>> added_pairs;
      for (size_t index = 0; index < 60; ++index) {
         const double infeasibility = std::pow(10., 4. * uniform(generator) - 2.);
         const double objective = 1e-3 * (uniform(generator) - 0.5);
         added_pairs.emplace_back(infeasibility, objective);
         filter.add(infeasibility, objective);
      }
      for (size_t sample = 0; sample < 200; ++sample) {
         const double infeasibility = std::pow(10., 4. * uniform(generator) - 2.);
         const double objective = 1e-3 * (uniform(generator) - 0.5) - 1e-2 * uniform(generator) * entry_margin_gamma * infeasibility;
         const bool taboo = in_entry_margin_taboo_region(added_pairs, filter_beta, entry_margin_gamma, infeasibility, objective);
         const bool acceptable = filter.filter_acceptable(infeasibility, objective);
         if (taboo && acceptable) {
            ++number_wrongly_accepted;
         }
         if (!taboo && !acceptable) {
            ++number_wrongly_rejected;
         }
      }
   }
   EXPECT_EQ(number_wrongly_accepted, 0u);
   EXPECT_EQ(number_wrongly_rejected, 0u);
}
