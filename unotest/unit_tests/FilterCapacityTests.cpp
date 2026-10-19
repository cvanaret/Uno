// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Fixed and unlimited filter capacities.
// Regression test: when the filter is full, dropping an entry must not make acceptable the points it rejected.
// Bug: the upper bound was set to beta * max(ub, h_last) = beta * ub, so the dropped entry's region was lost.

#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/switching_methods/filter_methods/Filter.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"

using namespace uno;

namespace {
   Options filter_options(size_t capacity) {
      Options options;
      DefaultOptions::load(options);
      options.set_integer("filter_capacity", static_cast<uno_int>(capacity));
      return options; // beta = 0.999, gamma = 0.001
   }

   bool is_acceptable(const Filter& filter, double infeasibility, double objective) {
      return filter.acceptable_wrt_infeasibility_upper_bound(infeasibility) &&
         filter.filter_acceptable(infeasibility, objective);
   }
} // namespace

// the entry with the largest infeasibility is dropped to make room for a new one
TEST(FilterCapacity, DroppedEntryRegionStaysRejected) {
   Filter filter(filter_options(2));
   filter.set_infeasibility_upper_bound(100.);
   filter.add(1., 10.);
   filter.add(2., 5.);
   ASSERT_FALSE(is_acceptable(filter, 3., 6.)); // dominated by (2, 5)

   filter.add(0.5, 20.); // full: (2, 5) is dropped
   EXPECT_FALSE(is_acceptable(filter, 3., 6.));
}

// the new pair has the largest infeasibility: the entries already in the filter must be kept
TEST(FilterCapacity, NewPairWithLargestInfeasibility) {
   Filter filter(filter_options(2));
   filter.set_infeasibility_upper_bound(100.);
   filter.add(1., 10.);
   filter.add(2., 5.);
   ASSERT_FALSE(is_acceptable(filter, 2.5, 6.)); // dominated by (2, 5)

   filter.add(3., 1.); // full
   EXPECT_FALSE(is_acceptable(filter, 2.5, 6.)); // still dominated by (2, 5)
   EXPECT_FALSE(is_acceptable(filter, 4., 2.)); // dominated by (3, 1)
   EXPECT_TRUE(is_acceptable(filter, 0.5, 0.)); // non-dominated pairs are still accepted
}

// below capacity, the upper bound is untouched
TEST(FilterCapacity, UpperBoundUnchangedWhenNotFull) {
   Filter filter(filter_options(3));
   filter.set_infeasibility_upper_bound(100.);
   filter.add(1., 10.);
   filter.add(2., 5.);
   EXPECT_TRUE(is_acceptable(filter, 50., -1e6));
}

// reset empties the filter and restores the upper bound lowered by the dropped entries
TEST(FilterCapacity, ResetRestoresUpperBound) {
   Filter filter(filter_options(2));
   filter.set_infeasibility_upper_bound(100.);
   filter.add(1., 10.);
   filter.add(2., 5.);
   filter.add(0.5, 20.); // full: (2, 5) is dropped and the upper bound is lowered to 2
   ASSERT_FALSE(is_acceptable(filter, 50., -1e6));

   filter.reset();
   EXPECT_TRUE(is_acceptable(filter, 50., -1e6));
   EXPECT_FALSE(is_acceptable(filter, 101., -1e6));
}

// capacity 0: the filter grows, no entry is dropped and the upper bound is never lowered
TEST(FilterCapacity, UnlimitedCapacity) {
   constexpr size_t number_pairs = 200;
   Filter unlimited_filter(filter_options(0));
   Filter limited_filter(filter_options(2));
   for (Filter* filter: {&unlimited_filter, &limited_filter}) {
      filter->set_infeasibility_upper_bound(1000.);
      for (size_t index = 1; index <= number_pairs; ++index) {
         filter->add(static_cast<double>(index), static_cast<double>(number_pairs - index));
      }
   }
   // every entry is kept: points slightly dominated by any entry are rejected
   for (size_t index = 1; index < number_pairs; ++index) {
      const double infeasibility = static_cast<double>(index) + 0.5;
      const double objective = static_cast<double>(number_pairs - index) + 0.1;
      EXPECT_FALSE(is_acceptable(unlimited_filter, infeasibility, objective));
   }
   EXPECT_TRUE(is_acceptable(unlimited_filter, 500., -1e6)); // upper bound untouched
   EXPECT_FALSE(is_acceptable(limited_filter, 500., -1e6)); // upper bound lowered by the dropped entries
   EXPECT_TRUE(is_acceptable(unlimited_filter, 0.5, 1e6)); // left of all entries
}
