// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <algorithm>
#include "Filter.hpp"
#include "tools/Logger.hpp"
#include "tools/Range.hpp"

Filter::Filter(const Options& options) :
      capacity(options.get_unsigned_int("filter_capacity")),
      infeasibility(this->capacity),
      optimality(this->capacity),
      parameters({
         options.get_double("filter_beta"),
         options.get_double("filter_gamma")
      }) {
   this->reset();
}

void Filter::reset() {
   this->upper_bound = INF<double>;
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
      return this->parameters.beta * this->upper_bound;
   }
}

void Filter::left_shift(size_t start, size_t shift_size) {
   for (size_t position: Range(start, this->number_entries - shift_size)) {
      this->infeasibility[position] = this->infeasibility[position + shift_size];
      this->optimality[position] = this->optimality[position + shift_size];
   }
}

void Filter::right_shift(size_t start, size_t shift_size) {
   for (size_t position: Range<BACKWARD>(this->number_entries, start)) {
      this->infeasibility[position] = this->infeasibility[position - shift_size];
      this->optimality[position] = this->optimality[position - shift_size];
   }
}

//  add (infeasibility_measure, optimality_measure) to the filter
void Filter::add(double infeasibility_measure, double optimality_measure) {
   // remove dominated filter entries
   // find position in filter without margin
   size_t start_position = 0;
   while (start_position < this->number_entries && this->infeasibility[start_position] < infeasibility_measure) {
      start_position++;
   }

   // find redundant entries starting from position
   size_t end_position = start_position;
   while (end_position < this->number_entries && optimality_measure <= this->optimality[end_position]) {
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
      this->upper_bound = this->parameters.beta * std::max(this->upper_bound, this->infeasibility[this->number_entries - 1]);
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
   this->optimality[start_position] = optimality_measure;
   this->number_entries++;
}

bool Filter::acceptable_wrt_upper_bound(double infeasibility_measure) const {
   return (infeasibility_measure < this->parameters.beta * this->upper_bound);
}

// return true if (infeasibility_measure, optimality_measure) acceptable, false otherwise
bool Filter::acceptable(double infeasibility_measure, double optimality_measure) {
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
   // until here, the optimality measure was not required
   else if (optimality_measure <= this->optimality[position - 1] - this->parameters.gamma * infeasibility_measure) {
      return true; // point acceptable
   }
   DEBUG << "Rejected because of filter domination\n";
   return false;
}

//! check acceptability wrt current point
bool Filter::acceptable_wrt_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
      double trial_optimality_measure) {
   return (trial_optimality_measure <= current_optimality_measure - this->parameters.gamma * trial_infeasibility_measure) ||
          (trial_infeasibility_measure < this->parameters.beta * current_infeasibility_measure);
}

double Filter::compute_actual_reduction(double current_optimality_measure, double /*current_infeasibility_measure*/, double trial_optimality_measure) {
   return current_optimality_measure - trial_optimality_measure;
}

//! print: print the content of the filter
std::ostream& operator<<(std::ostream& stream, Filter& filter) {
   stream << "************\n";
   stream << "  Current filter (infeasibility, optimality):\n";
   for (size_t position: Range(filter.number_entries)) {
      stream << "\t" << filter.infeasibility[position] << "\t" << filter.optimality[position] << '\n';
   }
   stream << "************\n";
   return stream;
}