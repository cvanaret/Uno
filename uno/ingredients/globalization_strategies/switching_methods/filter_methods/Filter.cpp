// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <stdexcept>
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
         infeasibility(this->capacity),
         objective(this->capacity),
         parameters({options.get_double("filter_beta"), options.get_double("filter_gamma")}) {
      if (this->capacity == 0) {
         throw std::runtime_error("The filter has capacity 0");
      }
   }

   void Filter::reset() {
      this->number_entries = 0;
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
   }

   bool Filter::acceptable_wrt_infeasibility_upper_bound(double trial_infeasibility) const {
      return this->infeasibility_sufficient_reduction(this->infeasibility_upper_bound, trial_infeasibility);
   }

   // return true if (infeasibility, objective) acceptable, false otherwise
   // note: the infeasibility upper bound must be tested separately (see acceptable_wrt_infeasibility_upper_bound)
   bool Filter::filter_acceptable(double trial_infeasibility, double trial_objective) const {
      // TODO: use binary search (use some form of https://en.cppreference.com/w/cpp/algorithm/binary_search.html)
      size_t position = 0;
      while (position < this->number_entries && !this->infeasibility_sufficient_reduction(this->infeasibility[position], trial_infeasibility)) {
         ++position;
      }

      // check acceptability
      if (position == 0) {
         return true; // acceptable as left-most entry
      }
      // until here, the objective measure was not evaluated
      else if (this->objective_sufficient_reduction(this->objective[position - 1], trial_objective, trial_infeasibility)) {
         return true; // point acceptable
      }
      DEBUG << "Rejected because of filter domination\n";
      return false;
   }

   //! check acceptability wrt current point
   bool Filter::acceptable_wrt_current_iterate(double current_infeasibility, double current_objective, double trial_infeasibility,
         double trial_objective) const {
      return this->infeasibility_sufficient_reduction(current_infeasibility, trial_infeasibility) ||
         this->objective_sufficient_reduction(current_objective, trial_objective, trial_infeasibility);
   }

   bool Filter::infeasibility_sufficient_reduction(double current_infeasibility, double trial_infeasibility) const {
      return (trial_infeasibility <= this->parameters.beta * current_infeasibility);
   }

   double Filter::compute_actual_objective_reduction(double current_objective, double /*current_infeasibility*/, double trial_objective) {
      return current_objective - trial_objective;
   }

   static std::string to_string(double number) {
      std::ostringstream stream;
      stream << std::defaultfloat << std::setprecision(7) << number;
      return stream.str();
   }

   // add (infeasibility, objective) to the filter
   // add (infeasibility, objective) to the filter
// invariant: infeasibility strictly increasing, objective strictly decreasing (raw values, no margins)
void Filter::add(double current_infeasibility, double current_objective) {
   // first entry whose infeasibility is >= the new one (raw comparison, same as the removal below)
   size_t position = 0;
   while (position < this->number_entries && this->infeasibility[position] < current_infeasibility) {
      ++position;
   }

   // if an existing entry dominates the new pair, the filter is unchanged:
   // - entries before `position` have smaller infeasibility; the one with the smallest objective is position-1
   // - the entry at `position` may have equal infeasibility
   if (0 < position && this->objective[position - 1] <= current_objective) {
      return;
   }
   if (position < this->number_entries && this->infeasibility[position] == current_infeasibility &&
         this->objective[position] <= current_objective) {
      return;
   }

   // remove the entries dominated by the new pair: from `position` on, infeasibility >= new one, and the
   // objectives are decreasing, so the dominated entries form a contiguous block
   size_t end_position = position;
   while (end_position < this->number_entries && current_objective <= this->objective[end_position]) {
      ++end_position;
   }
   const size_t number_dominated_entries = end_position - position;
   if (0 < number_dominated_entries) {
      this->left_shift(position, number_dominated_entries);
      this->number_entries -= number_dominated_entries;
   }

   // make room if the filter is full: drop the entry with the largest infeasibility and tighten the upper bound
   if (this->number_entries >= this->capacity) {
      const double largest_filter_infeasibility = std::max(this->infeasibility_upper_bound,
         this->infeasibility[this->number_entries - 1]);
      this->set_infeasibility_upper_bound(this->parameters.beta * largest_filter_infeasibility);
      --this->number_entries;
      position = std::min(position, this->number_entries);
   }

   // insert at `position` (no recomputation, no margin)
   if (position < this->number_entries) {
      this->right_shift(position, 1);
   }
   this->infeasibility[position] = current_infeasibility;
   this->objective[position] = current_objective;
   ++this->number_entries;
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
      for (size_t position: Range(filter.number_entries)) {
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

   bool Filter::is_empty() const {
      return (this->number_entries == 0);
   }

   bool Filter::objective_sufficient_reduction(double current_objective, double trial_objective, double trial_infeasibility) const {
      return (trial_objective <= current_objective - this->parameters.gamma * trial_infeasibility);
   }

   void Filter::left_shift(size_t start, size_t shift_size) {
      for (size_t position: Range(start, this->number_entries - shift_size)) {
         this->infeasibility[position] = this->infeasibility[position + shift_size];
         this->objective[position] = this->objective[position + shift_size];
      }
   }

   void Filter::right_shift(size_t start, size_t shift_size) {
      for (size_t position: BackwardRange(this->number_entries, start)) {
         this->infeasibility[position] = this->infeasibility[position - shift_size];
         this->objective[position] = this->objective[position - shift_size];
      }
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
      for (size_t index = 1; index < this->number_entries; ++index) {
         if (!(this->infeasibility[index - 1] < this->infeasibility[index] &&
               this->objective[index - 1] > this->objective[index])) {
            return false;
               }
      }
      return true;
   }
} // namespace