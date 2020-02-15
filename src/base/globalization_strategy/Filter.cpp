#include <iostream>
#include <cmath>
#include "Filter.hpp"

Filter::Filter(FilterConstants& constants) : constants(constants) {
    this->reset();
}

Filter::~Filter() {
}

void Filter::reset() {
    /* initialize the maximum filter size (not critical) */
    this->max_size = 50;
    this->upper_bound = INFINITY;
    this->entries.clear();
    return;
}

/*  add (infeasibility_measure, optimality_measure) to the filter */
void Filter::add(double infeasibility_measure, double optimality_measure) {
    /* remove dominated filter entries */
    std::list<FilterEntry>::iterator entry = this->entries.begin();
    while (entry != this->entries.end()) {
        if (infeasibility_measure < entry->infeasibility_measure && optimality_measure <= entry->optimality_measure) {
            entry = this->entries.erase(entry);
        } else {
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
    std::list<FilterEntry>::iterator position = this->entries.begin();
    while (position != this->entries.end() && infeasibility_measure >= this->constants.Beta * position->infeasibility_measure) {
        position++;
    }
    FilterEntry new_entry{infeasibility_measure, optimality_measure};
    this->entries.insert(position, new_entry);

    return;
}

/* query: return true if (infeasibility_measure, optimality_measure) acceptable, false otherwise */
bool Filter::query(double infeasibility_measure, double optimality_measure) {
    /* check upper bound first */
    if (this->constants.Beta * this->upper_bound <= infeasibility_measure) {
        return false;
    }

    std::list<FilterEntry>::iterator position = this->entries.begin();
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

//! query_current_iterate: check acceptable wrt current point 

bool Filter::query_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure, double trial_optimality_measure) {
    return !((trial_optimality_measure > current_optimality_measure - this->constants.Gamma * trial_infeasibility_measure) && (trial_infeasibility_measure > this->constants.Beta * current_infeasibility_measure));
}

double Filter::compute_actual_reduction(double current_objective, double current_residual, double trial_objective) {
    return current_objective - trial_objective;
}

//! print: print the content of the filter
std::ostream& operator<<(std::ostream &stream, Filter& filter) {
    stream << "************\n";
    stream << "  Current filter (constraint residual, objective):\n";
    for (FilterEntry const& entry : filter.entries) {
        stream << "\t" << entry.infeasibility_measure << "\t" << entry.optimality_measure << "\n";
    }
    stream << "************\n";
    return stream;
}

/* NonmonotoneFilter class*/

NonmonotoneFilter::NonmonotoneFilter(FilterConstants& constants, int number_dominated_entries) : Filter(constants), number_dominated_entries(number_dominated_entries) {
}

//!  add (infeasibility_measure,optimality_measure) to the filter

void NonmonotoneFilter::add(double infeasibility_measure, double optimality_measure) {
//    int new_size = this->entries.size();
//    int dominated_entries;
//    /* find entries in filter that are dominated by "number_dominated_entries" other entries */
//    for (int i = 0; i < this->size; i++) {
//        /* check whether ith entry dominated by (infeasibility_measure, optimality_measure) */
//        if ((this->optimality_measures[i] > optimality_measure) && (this->infeasibility_measures[i] > infeasibility_measure)) {
//            dominated_entries = 1;
//        }
//        else {
//            dominated_entries = 0;
//        }
//
//        /* find other filter entries that dominate ith entry */
//        for (int j = 0; j < this->size; j++) {
//            if ((this->optimality_measures[i] > this->optimality_measures[j]) && (this->infeasibility_measures[i] > this->infeasibility_measures[j])) {
//                dominated_entries++;
//            }
//        }
//        if (dominated_entries > this->number_dominated_entries) {
//            /* remove this entry */
//            this->shift_left(i, 1);
//            new_size--;
//        }
//    }
//
//    /* check sufficient space available */
//    if (new_size >= this->max_size) {
//        /* create space in filter: remove entry 1 (oldest entry) */
//        this->shift_left(1, 1);
//        new_size--;
//    }
//
//    /* add new entry to filter in position new_size */
//    this->infeasibility_measures[new_size] = infeasibility_measure;
//    this->optimality_measures[new_size] = optimality_measure;
//    this->size = new_size + 1;

    return;
}

//! query: check if (infeasibility_measure,optimality_measure) acceptable; return 1 if yes, 0 if not

bool NonmonotoneFilter::query(double infeasibility_measure, double optimality_measure) {
    /* check upper bound first */
    if (infeasibility_measure >= this->constants.Beta * this->upper_bound) {
        return false;
    }

    /* check acceptability by counting how many entries dominate */
    int dominated_entries = 0;
    for (FilterEntry const& entry : this->entries) {
        if (((optimality_measure > entry.optimality_measure - this->constants.Gamma * infeasibility_measure) && (infeasibility_measure >= this->constants.Beta * entry.infeasibility_measure)) ||
                ((optimality_measure >= entry.optimality_measure - this->constants.Gamma * infeasibility_measure) && (infeasibility_measure > this->constants.Beta * entry.infeasibility_measure))) {
            dominated_entries++;
        }
    }

    return (dominated_entries <= this->number_dominated_entries); // point acceptable (dominated by <= this->number_dominated_entries)
}

//! query_current_iterate: check acceptable wrt current point 

bool NonmonotoneFilter::query_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure, double trial_optimality_measure) {
    int dominated_entries;

    /* check acceptability wrt current point (non-monotone) */
    if ((trial_optimality_measure > current_optimality_measure - this->constants.Gamma * trial_infeasibility_measure) && (trial_infeasibility_measure > this->constants.Beta * current_infeasibility_measure)) {
        dominated_entries = 1;
    }
    else {
        dominated_entries = 0;
    }

    for (FilterEntry const& entry : this->entries) {
        if (((trial_optimality_measure > entry.optimality_measure - this->constants.Gamma * trial_infeasibility_measure) &&
                (trial_infeasibility_measure >= this->constants.Beta * entry.infeasibility_measure)) ||
                ((trial_optimality_measure >= entry.optimality_measure - this->constants.Gamma * trial_infeasibility_measure) &&
                (trial_infeasibility_measure > this->constants.Beta * entry.infeasibility_measure))) {
            dominated_entries++;
        }
    }
    return (dominated_entries <= this->number_dominated_entries); // point acceptable (dominated by <= this->number_dominated_entries)
}

/* check NON-MONOTONE sufficient reduction condition */
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
