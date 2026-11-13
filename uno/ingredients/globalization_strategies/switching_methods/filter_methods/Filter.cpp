// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <cstddef>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include "Filter.hpp"
#include "options/Options.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"
#include "tools/Symbols.hpp"

namespace uno {
   Filter::Filter(const Options& options) :
         capacity(options.get_unsigned_int("filter_capacity")),
         parameters({options.get_double("filter_beta"), options.get_double("filter_gamma")}),
         margin_uses_entry_infeasibility(options.get_bool("filter_margin_uses_entry_infeasibility")) {
      // fixed capacity: allocate once (one more slot to insert before dropping)
      if (0 < this->capacity) {
         this->infeasibility.reserve(this->capacity + 1);
         this->objective.reserve(this->capacity + 1);
      }
   }

   void Filter::reset() {
      this->infeasibility.clear();
      this->objective.clear();
      this->infeasibility_upper_bound = this->initial_infeasibility_upper_bound;
   }

   double Filter::get_smallest_infeasibility() const {
      if (!this->is_empty()) {
         // left-most entry has the lowest infeasibility
         return this->infeasibility[0];
      }
      else { // empty filter
         return this->infeasibility_upper_bound;
      }
   }

   void Filter::set_infeasibility_upper_bound(double new_upper_bound) {
      this->infeasibility_upper_bound = new_upper_bound;
      this->initial_infeasibility_upper_bound = new_upper_bound;
   }

   double Filter::get_infeasibility_upper_bound() const {
      return this->infeasibility_upper_bound;
   }

   bool Filter::acceptable_wrt_infeasibility_upper_bound(double trial_infeasibility) const {
      return this->infeasibility_sufficient_reduction(this->infeasibility_upper_bound, trial_infeasibility);
   }

   // return true if (infeasibility, objective) acceptable, false otherwise
   // note: the infeasibility upper bound must be tested separately (see acceptable_wrt_infeasibility_upper_bound)
   bool Filter::filter_acceptable(double trial_infeasibility, double trial_objective) const {
      // binary search: the infeasibility increases, so the entries that do not satisfy the infeasibility sufficient
      // reduction condition form a prefix. position: first entry that satisfies it
      const auto first_entry = this->infeasibility.cbegin();
      const size_t position = static_cast<size_t>(std::partition_point(first_entry, this->infeasibility.cend(),
         [&](double entry_infeasibility) {
            return !this->infeasibility_sufficient_reduction(entry_infeasibility, trial_infeasibility);
         }) - first_entry);

      // check acceptability
      if (position == 0) {
         return true; // acceptable as left-most entry
      }
      // until here, the objective measure was not evaluated
      if (this->objective_sufficient_reduction(this->objective[position - 1], trial_objective,
            this->margin_uses_entry_infeasibility ? this->infeasibility[position - 1] : trial_infeasibility)) {
         return true; // point acceptable
      }
      DEBUG << "Rejected because of filter domination\n";
      return false;
   }

   //! check acceptability wrt current point
   bool Filter::acceptable_wrt_current_iterate(double current_infeasibility, double current_objective, double trial_infeasibility,
         double trial_objective) const {
      return this->infeasibility_sufficient_reduction(current_infeasibility, trial_infeasibility) ||
         this->objective_sufficient_reduction(current_objective, trial_objective,
            this->margin_uses_entry_infeasibility ? current_infeasibility : trial_infeasibility);
   }

   bool Filter::infeasibility_sufficient_reduction(double current_infeasibility, double trial_infeasibility) const {
      return (trial_infeasibility <= this->parameters.beta * current_infeasibility);
   }

   static std::string to_string(double number) {
      std::ostringstream stream;
      stream << std::defaultfloat << std::setprecision(7) << number;
      return stream.str();
   }

   // add (infeasibility, objective) to the filter
   // invariant: infeasibility strictly increasing, objective strictly decreasing (raw values, no margins)
   void Filter::add(double current_infeasibility, double current_objective) {
      assert(this->infeasibility.size() == this->objective.size());
      // first entry whose infeasibility is >= the new one (raw comparison, same as the removal below)
      const auto infeasibility_begin = this->infeasibility.begin();
      const auto objective_begin = this->objective.begin();
      const size_t position = static_cast<size_t>(std::lower_bound(infeasibility_begin, this->infeasibility.end(),
         current_infeasibility) - infeasibility_begin);

      // if an existing entry dominates the new pair, the filter is unchanged:
      // - entries before `position` have smaller infeasibility; the one with the smallest objective is position-1
      // - the entry at `position` may have equal infeasibility
      if (0 < position && this->objective[position - 1] <= current_objective) {
         return;
      }
      if (position < this->number_entries() && this->infeasibility[position] == current_infeasibility &&
            this->objective[position] <= current_objective) {
         return;
      }

      // remove the entries dominated by the new pair: from `position` on, infeasibility >= new one, and the
      // objectives are decreasing, so the dominated entries form a contiguous block [position, end_position)
      const auto start = static_cast<std::ptrdiff_t>(position);
      const auto end = std::partition_point(objective_begin + start, this->objective.end(),
         [&](double entry_objective) { return current_objective <= entry_objective; }) - objective_begin;
      this->infeasibility.erase(infeasibility_begin + start, infeasibility_begin + end);
      this->objective.erase(objective_begin + start, objective_begin + end);

      // insert at `position` (no recomputation, no margin)
      this->infeasibility.insert(this->infeasibility.begin() + start, current_infeasibility);
      this->objective.insert(this->objective.begin() + start, current_objective);

      // fixed capacity exceeded: drop the pair (existing or new) with the largest infeasibility h_max and lower the upper
      // bound to h_max. Since the bound test is h <= beta*ub, every point with h > beta*h_max is rejected, which
      // includes the whole region dominated by the dropped pair
      if (0 < this->capacity && this->capacity < this->infeasibility.size()) {
         this->infeasibility_upper_bound = std::min(this->infeasibility_upper_bound, this->infeasibility.back());
         this->infeasibility.pop_back();
         this->objective.pop_back();
      }
      assert(this->is_sorted());
   }

   // print the content of the filter
   std::ostream& operator<<(std::ostream& stream, const Filter& filter) {
      std::string table;
      // header
      table.append(symbols::top_left_corner);
      for ([[maybe_unused]] size_t _: Range(Filter::fixed_length_column1 + 1)) {
         table.append(symbols::hyphen);
      }
      table.append(symbols::top_tee);
      for ([[maybe_unused]] size_t _: Range(Filter::fixed_length_column2 + 1)) {
         table.append(symbols::hyphen);
      }
      table.append(symbols::top_right_corner);
      table.push_back('\n');

      table.append(symbols::pipe);
      table.append(" infeasibility ");
      table.append(symbols::pipe);
      table.append(" objective ");
      table.append(symbols::pipe);
      table.push_back('\n');

      table.append(symbols::left_tee);
      for ([[maybe_unused]] size_t _: Range(Filter::fixed_length_column1 + 1)) {
         table.append(symbols::hyphen);
      }
      table.append(symbols::cross);
      for ([[maybe_unused]] size_t _: Range(Filter::fixed_length_column2 + 1)) {
         table.append(symbols::hyphen);
      }
      table.append(symbols::right_tee);
      table.push_back('\n');

      // values
      for (size_t position: Range(filter.number_entries())) {
         Filter::print_line(table, to_string(filter.infeasibility[position]), to_string(filter.objective[position]));
      }
      // print upper bound
      Filter::print_line(table, to_string(filter.infeasibility_upper_bound), "-");

      // footer
      table.append(symbols::bottom_left_corner);
      for ([[maybe_unused]] size_t _: Range(Filter::fixed_length_column1 + 1)) {
         table.append(symbols::hyphen);
      }
      table.append(symbols::bottom_tee);
      for ([[maybe_unused]] size_t _: Range(Filter::fixed_length_column2 + 1)) {
         table.append(symbols::hyphen);
      }
      table.append(symbols::bottom_right_corner);
      table.push_back('\n');

      stream << table;
      return stream;
   }

   // protected member functions

   size_t Filter::number_entries() const {
      return this->infeasibility.size();
   }

   bool Filter::is_empty() const {
      return this->infeasibility.empty();
   }


   bool Filter::objective_sufficient_reduction(double current_objective, double trial_objective, double trial_infeasibility) const {
      return (trial_objective <= current_objective - this->parameters.gamma * trial_infeasibility);
   }

   void Filter::print_line(std::string& table, const std::string& infeasibility, const std::string& objective) {
      // compute lengths of columns
      const size_t infeasibility_length = infeasibility.size();
      const size_t number_infeasibility_spaces = (infeasibility_length < fixed_length_column1) ? fixed_length_column1 - infeasibility_length : 0;
      const size_t objective_length = objective.size();
      const size_t number_objective_spaces = (objective_length < fixed_length_column2) ? fixed_length_column2 - objective_length : 0;

      // print line
      table.append(symbols::pipe);
      table.push_back(' ');
      table.append(infeasibility);
      for ([[maybe_unused]] size_t k: Range(number_infeasibility_spaces)) {
         table.push_back(' ');
      }
      table.append(symbols::pipe);
      table.push_back(' ');
      table.append(objective);
      for ([[maybe_unused]] size_t k: Range(number_objective_spaces)) {
         table.push_back(' ');
      }
      table.append(symbols::pipe);
      table.push_back('\n');
   }

   bool Filter::is_sorted() const {
      for (size_t index = 1; index < this->number_entries(); ++index) {
         if (!(this->infeasibility[index - 1] < this->infeasibility[index] &&
               this->objective[index - 1] > this->objective[index])) {
            return false;
               }
      }
      return true;
   }
} // namespace