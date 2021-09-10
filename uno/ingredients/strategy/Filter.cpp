#include <iostream>
#include <cmath>
#include "Filter.hpp"

Filter::Filter(FilterConstants& constants) : constants(constants) {
   this->reset();
}

void Filter::reset() {
   /* initialize the maximum filter size (not critical) */
   this->upper_bound = INFINITY;
   this->entries.clear();
}

/*  add (infeasibility_measure, optimality_measure) to the filter */
void Filter::add(double infeasibility_measure, double optimality_measure) {
   /* remove dominated filter entries */
   auto entry = this->entries.begin();
   while (entry != this->entries.end()) {
      if (infeasibility_measure < entry->infeasibility_measure && optimality_measure <= entry->optimality_measure) {
         entry = this->entries.erase(entry);
      }
      else {
         entry++;
      }
   }

   /* check sufficient space available for new entry (remove last entry, if not) */
   if (this->max_size <= this->entries.size()) {
      FilterEntry& last_element = this->entries.back();
      this->upper_bound = this->constants.Beta * std::max(this->upper_bound, last_element.infeasibility_measure);
      this->entries.pop_back();
   }

   /* add new entry to the filter */
   auto position = this->entries.begin();
   while (position != this->entries.end() && infeasibility_measure >= this->constants.Beta * position->infeasibility_measure) {
      position++;
   }
   FilterEntry new_entry{infeasibility_measure, optimality_measure};
   this->entries.insert(position, new_entry);
}

// filter must be nonempty

double Filter::eta_min() {
   FilterEntry& last_element = this->entries.back();
   return last_element.infeasibility_measure;
}

// filter must be nonempty

double Filter::omega_min() {
   FilterEntry& last_element = this->entries.back();
   return last_element.optimality_measure;
}

/* query: return true if (infeasibility_measure, optimality_measure) acceptable, false otherwise */
bool Filter::accept(double infeasibility_measure, double optimality_measure) {
   /* check upper bound first */
   if (this->constants.Beta * this->upper_bound <= infeasibility_measure) {
      return false;
   }

   auto position = this->entries.begin();
   while (position != this->entries.end() && infeasibility_measure >= this->constants.Beta * position->infeasibility_measure) {
      position++;
   }

   /* check acceptability */
   if (position == this->entries.begin()) {
      return true; // acceptable as left-most entry
   }
   else if (optimality_measure <= std::prev(position)->optimality_measure - this->constants.Gamma * infeasibility_measure) {
      return true; // point acceptable
   }
   else {
      return false; // point rejected
   }
}

//! improves_current_iterate: check acceptable wrt current point 

bool Filter::improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
      double trial_optimality_measure) {
   return (trial_optimality_measure <= current_optimality_measure - this->constants.Gamma * trial_infeasibility_measure) ||
          (trial_infeasibility_measure < this->constants.Beta * current_infeasibility_measure);
}

double Filter::compute_actual_reduction(double current_objective, double /*current_residual*/, double trial_objective) {
   return current_objective - trial_objective;
}

//! print: print the content of the filter

std::ostream& operator<<(std::ostream& stream, Filter& filter) {
   stream << "************\n";
   stream << "  Current filter (constraint residual, evaluate_objective):\n";
   for (FilterEntry const& entry: filter.entries) {
      stream << "\t" << entry.infeasibility_measure << "\t" << entry.optimality_measure << "\n";
   }
   stream << "************\n";
   return stream;
}

/* NonmonotoneFilter class */

NonmonotoneFilter::NonmonotoneFilter(FilterConstants& constants, int number_dominated_entries) :
      Filter(constants), number_dominated_entries(number_dominated_entries) {
}

//! add (infeasibility_measure, optimality_measure) to the filter
void NonmonotoneFilter::add(double infeasibility_measure, double optimality_measure) {
   int dominated_entries;
   /* find entries in filter that are dominated by "number_dominated_entries" other entries */
   auto entry = this->entries.begin();
   auto position = this->entries.end();
   while (entry != this->entries.end()) {
      /* check whether ith entry dominated by (infeasibility_measure, optimality_measure) */
      if ((optimality_measure < entry->optimality_measure) && (infeasibility_measure < entry->infeasibility_measure)) {
         dominated_entries = 1;
      }
      else {
         dominated_entries = 0;
      }

      /* find other filter entries that dominate ith entry */
      for (FilterEntry const& other_entry: this->entries) {
         if ((other_entry.optimality_measure < entry->optimality_measure) && (other_entry.infeasibility_measure < entry->infeasibility_measure)) {
            dominated_entries++;
         }
      }
      if (dominated_entries > this->number_dominated_entries) {
         /* remove this entry */
         entry = this->entries.erase(entry);
         position--;
      }
      else {
         entry++;
      }
   }

   /* check sufficient space available */
   if (this->max_size <= this->entries.size()) {
      /* create space in filter: remove entry 1 (oldest entry) */
      this->entries.pop_front();
   }

   /* add new entry to the filter */
   FilterEntry new_entry{infeasibility_measure, optimality_measure};
   this->entries.insert(position, new_entry);
}

//! accept: check if (infeasibility_measure, optimality_measure) acceptable
bool NonmonotoneFilter::accept(double infeasibility_measure, double optimality_measure) {
   /* check upper bound first */
   if (infeasibility_measure >= this->constants.Beta * this->upper_bound) {
      return false;
   }

   /* check acceptability by counting how many entries dominate */
   int dominated_entries = 0;
   for (FilterEntry const& entry: this->entries) {
      if ((optimality_measure > entry.optimality_measure - this->constants.Gamma * infeasibility_measure) &&
          (infeasibility_measure >= this->constants.Beta * entry.infeasibility_measure)) {
         dominated_entries++;
      }
      else if ((optimality_measure >= entry.optimality_measure - this->constants.Gamma * infeasibility_measure) &&
               (infeasibility_measure > this->constants.Beta * entry.infeasibility_measure)) {
         dominated_entries++;
      }
   }

   return (dominated_entries <= this->number_dominated_entries); // point acceptable (dominated by <= this->number_dominated_entries)
}

//! improves_current_iterate: check acceptable wrt current point
bool NonmonotoneFilter::improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure,
      double trial_infeasibility_measure, double trial_optimality_measure) {
   int dominated_entries;

   /* check acceptability wrt current point (non-monotone) */
   if ((trial_optimality_measure > current_optimality_measure - this->constants.Gamma * trial_infeasibility_measure) &&
       (trial_infeasibility_measure > this->constants.Beta * current_infeasibility_measure)) {
      dominated_entries = 1;
   }
   else {
      dominated_entries = 0;
   }

   for (FilterEntry const& entry: this->entries) {
      if ((trial_optimality_measure > entry.optimality_measure - this->constants.Gamma * trial_infeasibility_measure) &&
          (trial_infeasibility_measure >= this->constants.Beta * entry.infeasibility_measure)) {
         dominated_entries++;
      }
      else if ((trial_optimality_measure >= entry.optimality_measure - this->constants.Gamma * trial_infeasibility_measure) &&
               (trial_infeasibility_measure > this->constants.Beta * entry.infeasibility_measure)) {
         dominated_entries++;
      }
   }
   return (dominated_entries <= this->number_dominated_entries); // point acceptable (dominated by <= this->number_dominated_entries)
}

/* compute_actual_reduction: check nonmonotone sufficient reduction condition */
double NonmonotoneFilter::compute_actual_reduction(double current_objective, double current_residual, double trial_objective) {
   /* max penalty among most recent entries */
   double max_objective = current_objective;

   std::list<FilterEntry>::iterator position = this->entries.end();
   for (int i = 0; i < this->number_dominated_entries; i++) {
      double gamma;
      if (current_residual < position->infeasibility_measure) {
         gamma = 1. / this->constants.Gamma;
      }
      else {
         gamma = this->constants.Gamma;
      }
      double dash_objective = position->optimality_measure + (position->infeasibility_measure - current_residual) * gamma;
      max_objective = std::max(max_objective, dash_objective);
      position--;
   }
   /* non-monotone actual reduction */
   return max_objective - trial_objective;
}

/* FilterFactory class */

std::unique_ptr<Filter> FilterFactory::create(const Options& options) {
   double beta = stod(options.at("filter_Beta"));
   double gamma = stod(options.at("filter_Gamma"));
   FilterConstants filter_constants = {beta, gamma};
   std::string filter_type = options.at("strategy");

   if (filter_type == "filter") {
      return std::make_unique<Filter>(filter_constants);
   }
   else if (filter_type == "nonmonotone-filter") {
      int number_dominated_entries = stoi(options.at("nonmonotone_filter_number_dominated_entries"));
      return std::make_unique<NonmonotoneFilter>(filter_constants, number_dominated_entries);
   }
   else {
      throw std::invalid_argument("Filter type " + filter_type + " does not exist");
   }
}
