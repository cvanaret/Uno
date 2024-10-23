// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include "Filter.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"

namespace uno {
   Filter::Filter(const Options& options) :
         capacity(options.get_unsigned_int("filter_capacity")),
         infeasibility(this->capacity),
         objective(this->capacity),
         parameters({options.get_double("filter_beta"), options.get_double("filter_gamma")}) {
   }

   void Filter::reset() {
      this->number_entries = 0;
   }

   bool Filter::is_empty() const {
      return (this->number_entries == 0);
   }

   double Filter::get_smallest_infeasibility() const {
      if (not this->is_empty()) {
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

   //  add (infeasibility, objective) to the filter
   void Filter::add(double current_infeasibility, double current_objective) {
      // remove dominated filter entries
      // find position in filter without margin
      size_t start_position = 0;
      while (start_position < this->number_entries && this->infeasibility[start_position] < current_infeasibility) {
         start_position++;
      }

      // find redundant entries starting from position
      size_t end_position = start_position;
      while (end_position < this->number_entries && current_objective <= this->objective[end_position]) {
         end_position++;
      }

      // remove entries [position:end_position] from filter
      const size_t number_redundant_entries = end_position - start_position;
      if (0 < number_redundant_entries) {
         this->left_shift(start_position, number_redundant_entries);
         this->number_entries -= number_redundant_entries;
      }

      // check sufficient space available for new entry (remove last entry, if not)
      if (this->number_entries >= this->capacity) {
         const double largest_filter_infeasibility = std::max(this->infeasibility_upper_bound, this->infeasibility[this->number_entries - 1]);
         this->set_infeasibility_upper_bound(this->parameters.beta * largest_filter_infeasibility);
         // create space in filter: remove last entry
         this->number_entries--;
      }

      // add new entry to the filter at position
      start_position = 0;
      while (start_position < this->number_entries && not this->infeasibility_sufficient_reduction(this->infeasibility[start_position], current_infeasibility)) {
         start_position++;
      }
      // shift entries by one to right to make room for new entry
      if (start_position < this->number_entries) {
         this->right_shift(start_position, 1);
      }
      // add new entry to filter
      this->infeasibility[start_position] = current_infeasibility;
      this->objective[start_position] = current_objective;
      this->number_entries++;
   }

   bool Filter::acceptable_wrt_upper_bound(double trial_infeasibility) const {
      return this->infeasibility_sufficient_reduction(this->infeasibility_upper_bound, trial_infeasibility);
   }

   // return true if (infeasibility, objective) acceptable, false otherwise
   bool Filter::acceptable(double trial_infeasibility, double trial_objective) {
      // check upper bound first
      if (not this->acceptable_wrt_upper_bound(trial_infeasibility)) {
         DEBUG << "Rejected because of filter upper bound\n";
         return false;
      }

      // TODO: use binary search
      size_t position = 0;
      while (position < this->number_entries && not this->infeasibility_sufficient_reduction(this->infeasibility[position], trial_infeasibility)) {
         position++;
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

   double Filter::compute_actual_objective_reduction(double current_objective, double /*current_infeasibility*/, double trial_objective) {
      return current_objective - trial_objective;
   }

   bool Filter::infeasibility_sufficient_reduction(double current_infeasibility, double trial_infeasibility) const {
      return (trial_infeasibility < this->parameters.beta * current_infeasibility);
   }

   bool Filter::objective_sufficient_reduction(double current_objective, double trial_objective, double trial_infeasibility) const {
      return (trial_objective <= current_objective - this->parameters.gamma * trial_infeasibility);
   }

   std::string to_string(double number) {
      std::ostringstream stream;
      stream << std::defaultfloat << std::setprecision(7) << number;
      return stream.str();
   }

   void print_line(std::ostream& stream, const std::string& infeasibility, const std::string& objective) {
      const size_t fixed_length_column1 = 14;
      const size_t fixed_length_column2 = 10;

      // compute lengths of columns
      const size_t infeasibility_length = infeasibility.size();
      const size_t number_infeasibility_spaces = (infeasibility_length < fixed_length_column1) ? fixed_length_column1 - infeasibility_length : 0;
      const size_t objective_length = objective.size();
      const size_t number_objective_spaces = (objective_length < fixed_length_column2) ? fixed_length_column2 - objective_length : 0;

      // print line
      stream << "│ " << infeasibility;
      for ([[maybe_unused]] size_t k: Range(number_infeasibility_spaces)) {
         stream << ' ';
      }
      stream << "│ " << objective;
      for ([[maybe_unused]] size_t k: Range(number_objective_spaces)) {
         stream << ' ';
      }
      stream << "│\n";
   }

   // print the content of the filter
   std::ostream& operator<<(std::ostream& stream, Filter& filter) {
      stream << "┌───────────────┬───────────┐\n";
      stream << "│ infeasibility │ objective │\n";
      stream << "├───────────────┼───────────┤\n";
      for (size_t position: Range(filter.number_entries)) {
         print_line(stream, to_string(filter.infeasibility[position]), to_string(filter.objective[position]));
      }
      // print upper bound
      print_line(stream, to_string(filter.infeasibility_upper_bound), "-");
      stream << "└───────────────┴───────────┘\n";
      return stream;
   }
} // namespace
