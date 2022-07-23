// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <algorithm>
#include "Filter.hpp"
#include "tools/Logger.hpp"

Filter::Filter(const Options& options) :
      capacity(stoul(options.at("filter_capacity"))),
      infeasibility(this->capacity),
      optimality(this->capacity),
      constants({stod(options.at("filter_beta")),
                 stod(options.at("filter_gamma"))}) {
   this->reset();
}

void Filter::reset() {
   // initialize the maximum filter size (not critical)
   this->upper_bound = INF;
   this->number_entries = 0;
}

void Filter::left_shift(size_t start, size_t shift_size) {
   for (size_t position = start; position < this->number_entries - shift_size; position++) {
      this->infeasibility[position] = this->infeasibility[position + shift_size];
      this->optimality[position] = this->optimality[position + shift_size];
   }
}

void Filter::right_shift(size_t start, size_t shift_size) {
   for (size_t position = this->number_entries; position > start; position--) {
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
   if (this->number_entries > this->capacity - 1) {
      this->upper_bound = this->constants.beta * std::max(this->upper_bound, this->infeasibility[this->number_entries - 1]);
      // create space in filter: remove last entry
      this->number_entries--;
   }

   // add new entry to the filter at position
   start_position = 0;
   while (start_position < this->number_entries && infeasibility_measure >= this->constants.beta * this->infeasibility[start_position]) {
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

// query: return true if (infeasibility_measure, optimality_measure) acceptable, false otherwise
bool Filter::accept(double infeasibility_measure, double optimality_measure) {
   // check upper bound first
   if (this->constants.beta * this->upper_bound <= infeasibility_measure) {
      DEBUG << "Rejected because of filter upper bound\n";
      return false;
   }

   size_t position = 0;
   while (position < this->number_entries && infeasibility_measure >= this->constants.beta * this->infeasibility[position]) {
      position++;
   }

   // check acceptability
   if (position == 0) {
      return true; // acceptable as left-most entry
   }
   else if (optimality_measure <= this->optimality[position - 1] - this->constants.gamma * infeasibility_measure) {
      return true; // point acceptable
   }
   DEBUG << "Rejected because the optimality measure is not low enough\n";
   return false;
}

//! improves_current_iterate: check acceptable wrt current point
bool Filter::improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
      double trial_optimality_measure) {
   return (trial_optimality_measure <= current_optimality_measure - this->constants.gamma * trial_infeasibility_measure) ||
          (trial_infeasibility_measure < this->constants.beta * current_infeasibility_measure);
}

double Filter::compute_actual_reduction(double current_objective, double /*current_residual*/, double trial_objective) {
   return current_objective - trial_objective;
}

//! print: print the content of the filter
std::ostream& operator<<(std::ostream& stream, Filter& filter) {
   stream << "************\n";
   stream << "  Current filter (infeasibility, optimality):\n";
   for (size_t position = 0; position < filter.number_entries; position++) {
      stream << "\t" << filter.infeasibility[position] << "\t" << filter.optimality[position] << '\n';
   }
   stream << "************\n";
   return stream;
}

// NonmonotoneFilter class
NonmonotoneFilter::NonmonotoneFilter(const Options& options) :
      Filter(options), max_number_dominated_entries(stoul(options.at("nonmonotone_filter_number_dominated_entries"))) {
}

//! add (infeasibility_measure, optimality_measure) to the filter
void NonmonotoneFilter::add(double infeasibility_measure, double optimality_measure) {
   // find entries in filter that are dominated by M other entries
   for (size_t i = 0; i < this->number_entries; i++) {
      size_t number_dominated = 0;
      // check whether ith entry dominated by (infeasibility_measure,optimality_measure)
      if ((this->optimality[i] > optimality_measure) && (this->infeasibility[i] > infeasibility_measure)) {
         number_dominated = 1;
      }
      // find other filter entries that dominate ith entry
      for (size_t j = 0; j < this->number_entries; j++) {
         if ((this->optimality[i] > this->optimality[j]) && (this->infeasibility[i] > this->infeasibility[j])) {
            number_dominated++;
         }
      }
      if (number_dominated > this->max_number_dominated_entries) {
         // remove this entry
         this->left_shift(i, 1);
         this->number_entries--;
      }
   }

   // check sufficient space available
   if (this->number_entries > this->capacity - 1) {
      // create space in filter: remove entry 1 (oldest entry)
      this->left_shift(1, 1);
      this->number_entries--;
   }

   // add new entry to filter in position this->number_entries
   this->infeasibility[this->number_entries] = infeasibility_measure;
   this->optimality[this->number_entries] = optimality_measure;
   this->number_entries++;
}

size_t NonmonotoneFilter::compute_number_dominated_entries(double infeasibility_measure, double optimality_measure) {
   size_t number_dominated_entries = 0;
   for (size_t i = 0; i < this->number_entries; i++) {
      if ((optimality_measure > this->optimality[i] - this->constants.gamma * infeasibility_measure) &&
          (infeasibility_measure >= this->constants.beta * this->infeasibility[i])) {
         number_dominated_entries++;
      }
      else if ((optimality_measure >= this->optimality[i] - this->constants.gamma * infeasibility_measure) &&
               (infeasibility_measure > this->constants.beta * this->infeasibility[i])) {
         number_dominated_entries++;
      }
   }
   return number_dominated_entries;
}

//! accept: check if (infeasibility_measure, optimality_measure) acceptable
bool NonmonotoneFilter::accept(double infeasibility_measure, double optimality_measure) {
   // check upper bound first
   if (infeasibility_measure >= this->constants.beta * this->upper_bound) {
      return false;
   }

   // check acceptability ** by counting how many entries dominate **
   const size_t number_dominated_entries = this->compute_number_dominated_entries(infeasibility_measure, optimality_measure);
   return (number_dominated_entries <= this->max_number_dominated_entries);
}

//! improves_current_iterate: check acceptable wrt current point
bool NonmonotoneFilter::improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure,
      double trial_infeasibility_measure, double trial_optimality_measure) {
   // check acceptability wrt current point (non-monotone)
   size_t number_dominated_entries = this->compute_number_dominated_entries(trial_infeasibility_measure, trial_optimality_measure);

   if ((trial_optimality_measure > current_optimality_measure - this->constants.gamma * trial_infeasibility_measure) &&
       (trial_infeasibility_measure > this->constants.beta * current_infeasibility_measure)) {
      number_dominated_entries++;
   }
   return (number_dominated_entries <= this->max_number_dominated_entries);
}

// compute_actual_reduction: check nonmonotone sufficient reduction condition
double NonmonotoneFilter::compute_actual_reduction(double current_objective, double current_residual, double trial_objective) {
   // check NON-MONOTONE sufficient reduction condition
   // max penalty among most recent entries
   double max_objective = current_objective;
   for (size_t i = 0; i < this->max_number_dominated_entries; i++) {
      const double gamma = (current_residual < this->infeasibility[this->number_entries - i]) ? 1 / this->constants.gamma : this->constants.gamma;
      const double dash_objective = this->optimality[this->number_entries - i] + gamma * (this->infeasibility[this->number_entries - i] -
         current_residual);
      max_objective = std::max(max_objective, dash_objective);
   }
   // non-monotone actual reduction
   return max_objective - trial_objective;
}

// FilterFactory class
std::unique_ptr<Filter> FilterFactory::create(const Options& options) {
   std::string filter_type = options.at("strategy");
   if (filter_type == "filter") {
      return std::make_unique<Filter>(options);
   }
   else if (filter_type == "nonmonotone-filter") {
      return std::make_unique<NonmonotoneFilter>(options);
   }
   else {
      throw std::invalid_argument("Filter type " + filter_type + " does not exist");
   }
}
