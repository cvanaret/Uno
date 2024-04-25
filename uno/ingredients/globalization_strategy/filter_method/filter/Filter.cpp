// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "Filter.hpp"
#include "tools/Logger.hpp"
#include "tools/Range.hpp"

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
      // left-most entry has the lowest infeasibility. Relax it with the envelope coefficient
      return this->parameters.beta * this->infeasibility[0];
   }
   else { // filter empty
      return this->parameters.beta * this->infeasibility_upper_bound;
   }
}

double Filter::get_infeasibility_upper_bound() const {
   return this->infeasibility_upper_bound;
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
   for (size_t position: Range<BACKWARD>(this->number_entries, start)) {
      this->infeasibility[position] = this->infeasibility[position - shift_size];
      this->objective[position] = this->objective[position - shift_size];
   }
}

//  add (infeasibility_measure, objective_measure) to the filter
void Filter::add(double infeasibility_measure, double objective_measure) {
   // remove dominated filter entries
   // find position in filter without margin
   size_t start_position = 0;
   while (start_position < this->number_entries && this->infeasibility[start_position] < infeasibility_measure) {
      start_position++;
   }

   // find redundant entries starting from position
   size_t end_position = start_position;
   while (end_position < this->number_entries && objective_measure <= this->objective[end_position]) {
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
      this->infeasibility_upper_bound = this->parameters.beta * std::max(this->infeasibility_upper_bound, this->infeasibility[this->number_entries - 1]);
      // create space in filter: remove last entry
      this->number_entries--;
   }

   // add new entry to the filter at position
   start_position = 0;
   while (start_position < this->number_entries && infeasibility_measure >= this->parameters.beta * this->infeasibility[start_position]) {
      start_position++;
   }
   // shift entries by one to right to make room for new entry
   if (start_position < this->number_entries) {
      this->right_shift(start_position, 1);
   }
   // add new entry to filter
   this->infeasibility[start_position] = infeasibility_measure;
   this->objective[start_position] = objective_measure;
   this->number_entries++;
}

bool Filter::acceptable_wrt_upper_bound(double infeasibility_measure) const {
   return (infeasibility_measure < this->parameters.beta * this->infeasibility_upper_bound);
}

// return true if (infeasibility_measure, objective_measure) acceptable, false otherwise
bool Filter::acceptable(double infeasibility_measure, double objective_measure) {
   // check upper bound first
   if (not this->acceptable_wrt_upper_bound(infeasibility_measure)) {
      DEBUG << "Rejected because of filter upper bound\n";
      return false;
   }

   // TODO: use binary search
   size_t position = 0;
   while (position < this->number_entries && infeasibility_measure >= this->parameters.beta * this->infeasibility[position]) {
      position++;
   }

   // check acceptability
   if (position == 0) {
      return true; // acceptable as left-most entry
   }
   // until here, the objective measure was not evaluated
   else if (objective_measure <= this->objective[position - 1] - this->parameters.gamma * infeasibility_measure) {
      return true; // point acceptable
   }
   DEBUG << "Rejected because of filter domination\n";
   return false;
}

//! check acceptability wrt current point
bool Filter::acceptable_wrt_current_iterate(double current_infeasibility_measure, double current_objective_measure, double trial_infeasibility_measure,
      double trial_objective_measure) {
   return (trial_objective_measure <= current_objective_measure - this->parameters.gamma * trial_infeasibility_measure) ||
          (trial_infeasibility_measure < this->parameters.beta * current_infeasibility_measure);
}

double Filter::compute_actual_objective_reduction(double current_objective_measure, double /*current_infeasibility_measure*/, double trial_objective_measure) {
   return current_objective_measure - trial_objective_measure;
}

std::string to_string(double number) {
   std::ostringstream stream;
   stream << std::defaultfloat << std::setprecision(7) << number;
   return stream.str();
}

// print the content of the filter
std::ostream& operator<<(std::ostream& stream, Filter& filter) {
   const size_t fixed_length_column1 = 14;
   const size_t fixed_length_column2 = 10;

   stream << "┌───────────────┬───────────┐\n";
   stream << "│ infeasibility │ objective │\n";
   stream << "├───────────────┼───────────┤\n";
   for (size_t position: Range(filter.number_entries)) {
      // convert numbers to strings
      const std::string infeasibility_string = to_string(filter.infeasibility[position]);
      const std::string objective_string = to_string(filter.objective[position]);

      // compute lengths of columns
      const size_t infeasibility_length = infeasibility_string.size();
      const size_t number_infeasibility_spaces = (infeasibility_length < fixed_length_column1) ? fixed_length_column1 - infeasibility_length : 0;
      const size_t objective_length = objective_string.size();
      const size_t number_objective_spaces = (objective_length < fixed_length_column2) ? fixed_length_column2 - objective_length : 0;

      // print line
      stream << "│ " << infeasibility_string;
      for ([[maybe_unused]] size_t k: Range(number_infeasibility_spaces)) {
         stream << ' ';
      }
      stream << "│ " << objective_string;
      for ([[maybe_unused]] size_t k: Range(number_objective_spaces)) {
         stream << ' ';
      }
      stream << "│\n";
   }
   stream << "└───────────────┴───────────┘\n";
   std::cout << "Infeasibility upper bound: " << filter.infeasibility_upper_bound << '\n';
   return stream;
}