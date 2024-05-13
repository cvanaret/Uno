// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "NonmonotoneFilter.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"

NonmonotoneFilter::NonmonotoneFilter(const Options& options) :
      Filter(options), max_number_dominated_entries(options.get_unsigned_int("nonmonotone_filter_number_dominated_entries")) {
}

//! add (infeasibility_measure, objective_measure) to the filter
void NonmonotoneFilter::add(double infeasibility_measure, double objective_measure) {
   // find entries in filter that are dominated by M other entries
   for (size_t entry_index: Range(this->number_entries)) {
      size_t number_dominated = 0;
      // check whether ith entry dominated by (infeasibility_measure,objective_measure)
      if ((this->objective[entry_index] > objective_measure) && (this->infeasibility[entry_index] > infeasibility_measure)) {
         number_dominated = 1;
      }
      // find other filter entries that dominate ith entry
      for (size_t other_entry_index: Range(this->number_entries)) {
         if ((this->objective[entry_index] > this->objective[other_entry_index]) && (this->infeasibility[entry_index] > this->infeasibility[other_entry_index])) {
            number_dominated++;
         }
      }
      if (number_dominated > this->max_number_dominated_entries) {
         // remove this entry
         this->left_shift(entry_index, 1);
         this->number_entries--;
      }
   }

   // check sufficient space available
   if (this->number_entries >= this->capacity) {
      // create space in filter: remove entry 1 (oldest entry)
      this->left_shift(1, 1);
      this->number_entries--;
   }

   // add new entry to filter in position this->number_entries
   this->infeasibility[this->number_entries] = infeasibility_measure;
   this->objective[this->number_entries] = objective_measure;
   this->number_entries++;
}

size_t NonmonotoneFilter::compute_number_dominated_entries(double infeasibility_measure, double objective_measure) {
   size_t number_dominated_entries = 0;
   for (size_t entry_index: Range(this->number_entries)) {
      if ((objective_measure > this->objective[entry_index] - this->parameters.gamma * infeasibility_measure) &&
          (infeasibility_measure >= this->parameters.beta * this->infeasibility[entry_index])) {
         number_dominated_entries++;
      }
      else if ((objective_measure >= this->objective[entry_index] - this->parameters.gamma * infeasibility_measure) &&
               (infeasibility_measure > this->parameters.beta * this->infeasibility[entry_index])) {
         number_dominated_entries++;
      }
   }
   return number_dominated_entries;
}

//! accept: check if (infeasibility_measure, objective_measure) acceptable
bool NonmonotoneFilter::acceptable(double infeasibility_measure, double objective_measure) {
   // check upper bound first
   if (not this->acceptable_wrt_upper_bound(infeasibility_measure)) {
      DEBUG << "Rejected because of filter upper bound\n";
      return false;
   }

   // check acceptability ** by counting how many entries dominate **
   // until here, the objective measure was not evaluated
   const size_t number_dominated_entries = this->compute_number_dominated_entries(infeasibility_measure, objective_measure);
   return (number_dominated_entries <= this->max_number_dominated_entries);
}

//! check acceptability wrt current point
bool NonmonotoneFilter::acceptable_wrt_current_iterate(double current_infeasibility_measure, double current_objective_measure,
      double trial_infeasibility_measure, double trial_objective_measure) {
   // check acceptability wrt current point (non-monotone)
   size_t number_dominated_entries = this->compute_number_dominated_entries(trial_infeasibility_measure, trial_objective_measure);

   if ((trial_objective_measure > current_objective_measure - this->parameters.gamma * trial_infeasibility_measure) &&
       (trial_infeasibility_measure > this->parameters.beta * current_infeasibility_measure)) {
      number_dominated_entries++;
   }
   return (number_dominated_entries <= this->max_number_dominated_entries);
}

// compute_actual_reduction: check nonmonotone sufficient reduction condition
double NonmonotoneFilter::compute_actual_objective_reduction(double current_objective_measure, double current_infeasibility_measure, double trial_objective_measure) {
   // check NON-MONOTONE sufficient reduction condition
   // max penalty among most recent entries
   double max_objective = current_objective_measure;
   for (size_t entry_index: Range(this->max_number_dominated_entries)) {
      const size_t index = this->number_entries - entry_index;
      const double gamma = (current_infeasibility_measure < this->infeasibility[index]) ? 1 / this->parameters.gamma : this->parameters.gamma;
      const double dash_objective = this->objective[index] + gamma * (this->infeasibility[index] - current_infeasibility_measure);
      max_objective = std::max(max_objective, dash_objective);
   }
   // non-monotone actual reduction
   return max_objective - trial_objective_measure;
}