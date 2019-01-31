#include <iostream>
#include <cmath>
#include "Filter.hpp"

Filter::Filter(FilterConstants& constants): constants(constants) {
	this->reset();
}

Filter::~Filter() {
}

void Filter::reset() {
	/* initialize the maximum filter size (not critical) */
	this->max_size = 50;
	this->size = 0;
	this->upper_bound = INFINITY;
	//this->upper_bound = std::numeric_limits<double>::infinity();

	/* allocate storage for filter storage */
	//this->infeasibility_measures.reserve(this->max_size);
	//this->optimality_measures.reserve(this->max_size);
    this->infeasibility_measures = std::vector<double>(this->max_size);
	this->optimality_measures = std::vector<double>(this->max_size);
	return;
}

/*  add (infeasibility_measure, optimality_measure) to the filter */
void Filter::add(double infeasibility_measure, double optimality_measure) {
	/* remove any redundant entries from the filter */

	/* find position in filter WITHOUT margin */
	int start = 0;
	while (start < this->size && infeasibility_measure >= this->infeasibility_measures[start]) {
		start++;
	}

	/* find redundant entries starting from start */
	int end = start;
	while (end < this->size && optimality_measure <= this->optimality_measures[end]) {
		end++;
	}
	--end;

	int trial_size = this->size;
	/* remove entries [start:end] from filter */
	int length = end - start + 1;
	if (0 < length) {
		this->shift_left(start, length);
		trial_size -= length;
	}

	/* check sufficient space available for new entry (remove last entry, if not) */
	if (this->max_size <= this->size) {
		/* set upper bound */
		this->upper_bound = this->constants.Beta*std::max(this->upper_bound, this->infeasibility_measures[this->size-1]);
		/* create space in filter: remove entry this->size-1 */
		trial_size--;
	}

	/* add new entry to the filter */

	/*  find position of infeasibility_measure in filter (new entry will be index) */
	int index = 0;
	while (index < trial_size && infeasibility_measure >= this->constants.Beta*this->infeasibility_measures[index]) {
		index++;
	}

	/* shift entries by one to right to make room for new entry */
	if (index < trial_size) {
		int length = 1;
		this->size = trial_size; // filter size changes above
		this->shift_right(index, length);
	}

	/* add new entry to filter in position index */
	this->infeasibility_measures[index] = infeasibility_measure;
	this->optimality_measures[index] = optimality_measure;
	this->size = trial_size + 1;

	return;
}

/* query: return true if (infeasibility_measure,optimality_measure) acceptable, false otherwise */
bool Filter::query(double infeasibility_measure, double optimality_measure) {
	/* check upper bound first */
	if (this->constants.Beta*this->upper_bound <= infeasibility_measure) {
		return false;
	}

	/* find position of infeasibility_measure in filter (new entry will be index) */
	int index = 0;
	while (index < this->size && infeasibility_measure >= this->constants.Beta*this->infeasibility_measures[index]) {
		index++;
	}

	/* now check acceptability */
	if (index == 0)
		return true; // acceptable as left-most entry
	else if (optimality_measure <= this->optimality_measures[index-1] - this->constants.Gamma*infeasibility_measure) 
		return true; // point acceptable
	else 
		return false; // point rejected by entry [index-1]
}

//! query_current_iterate: check acceptable wrt current point 
bool Filter::query_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure, double trial_optimality_measure) {
	return !((trial_optimality_measure > current_optimality_measure - this->constants.Gamma*trial_infeasibility_measure) && (trial_infeasibility_measure > this->constants.Beta*current_infeasibility_measure));
}

//! shift_left: shift filter to left from start to end
void Filter::shift_left(int start, int length) {
	for (int i = start; i < this->size-length; i++) {
		this->infeasibility_measures[i] = this->infeasibility_measures[i+length];
		this->optimality_measures[i] = this->optimality_measures[i+length];
	}
	return;
}

//! shift_right: shift filter to right from end to start
void Filter::shift_right(int start, int length) {
	for (int i = this->size; i > start; i--) {
		this->infeasibility_measures[i] = this->infeasibility_measures[i-length];
		this->optimality_measures[i] = this->optimality_measures[i-length];
	}
	return;
}

double Filter::compute_actual_reduction(double current_objective, double current_residual, double trial_objective) {
	return current_objective - trial_objective;
}

//! print: print the complete filter
std::ostream& operator<< (std::ostream &stream, Filter& filter) {
	stream << "    Current filter (constraint residual, objective):\n";
	for (int i = 0; i < filter.size; i++) {
		//fprintf(stdout, "\t%16.8g %16.8g\n", this->infeasibility_measures[i], this->optimality_measures[i]);
		stream << "\t" << filter.infeasibility_measures[i] << "\t" << filter.optimality_measures[i] << "\n";
	}
	//fprintf(stdout, "    Upper bound = %16.8g\n", this->upper_bound);
	stream << "    Upper bound = " << filter.upper_bound << "\n";
	return stream;
}

/* NonmonotoneFilter class*/

NonmonotoneFilter::NonmonotoneFilter(FilterConstants& constants, int number_dominated_entries):
		Filter(constants), number_dominated_entries(number_dominated_entries) {
}

//!  add (infeasibility_measure,optimality_measure) to the filter
void NonmonotoneFilter::add(double infeasibility_measure, double optimality_measure) {
	int trial_size = this->size;
	int dominated_entries;
	/* find entries in filter that are dominated by this->number_dominated_entries other entries */
	for (int i = 0; i < this->size; i++) {
		/* check whether ith entry dominated by (infeasibility_measure,optimality_measure) */
		if ((this->optimality_measures[i] > optimality_measure) && (this->infeasibility_measures[i] > infeasibility_measure)) {
			dominated_entries = 1;
		}
		else {
			dominated_entries = 0;
		}
		
		/* find other filter entries that dominate ith entry */
		for (int j = 0; j < this->size; j++) {
			if ((this->optimality_measures[i] > this->optimality_measures[j]) && (this->infeasibility_measures[i] > this->infeasibility_measures[j]))  {
				dominated_entries++;
			}
		}
		if (dominated_entries > this->number_dominated_entries) {
			/* remove this entry */
			this->shift_left(i, 1);
			trial_size--;
		}
	}

	/* check sufficient space available */
	if (trial_size > this->max_size - 1) {
		/* create space in filter: remove entry 1 (oldest entry) */
		this->shift_left(1, 1);
		trial_size--;
	}

	/* add new entry to filter in position trial_size */
	this->infeasibility_measures[trial_size] = infeasibility_measure;
	this->optimality_measures[trial_size] = optimality_measure;
	this->size = trial_size + 1;

	return;
}

//! query: check if (infeasibility_measure,optimality_measure) acceptable; return 1 if yes, 0 if not
bool NonmonotoneFilter::query(double infeasibility_measure, double optimality_measure) {
	/* check upper bound first */
	if (infeasibility_measure >= this->constants.Beta*this->upper_bound) {
		return 0;
	}
	
	/* check acceptability by counting how many entries dominate */
	int dominated_entries = 0;
	for (int i = 0; i < this->size; i++) {
		if (((optimality_measure > this->optimality_measures[i] - this->constants.Gamma*infeasibility_measure) && (infeasibility_measure >= this->constants.Beta*this->infeasibility_measures[i])) ||
				((optimality_measure >= this->optimality_measures[i] - this->constants.Gamma*infeasibility_measure) && (infeasibility_measure > this->constants.Beta*this->infeasibility_measures[i]))) {
			dominated_entries++;
		}
	}
	
	return (dominated_entries <= this->number_dominated_entries); // point acceptable (dominated by <= this->number_dominated_entries)
}

//! query_current_iterate: check acceptable wrt current point 
bool NonmonotoneFilter::query_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure, double trial_optimality_measure) {
	int dominated_entries;
	
	/* check acceptability wrt current point (non-monotone) */
	if ((trial_optimality_measure > current_optimality_measure - this->constants.Gamma*trial_infeasibility_measure) && (trial_infeasibility_measure > this->constants.Beta*current_infeasibility_measure)) 
		dominated_entries = 1;
	else
		dominated_entries = 0;
		
	for (int i = 0; i < this->size; i++) {
		if (((trial_optimality_measure > this->optimality_measures[i] - this->constants.Gamma*trial_infeasibility_measure) &&
                (trial_infeasibility_measure >= this->constants.Beta*this->infeasibility_measures[i])) ||
                ((trial_optimality_measure >= this->optimality_measures[i] - this->constants.Gamma*trial_infeasibility_measure) &&
                (trial_infeasibility_measure > this->constants.Beta*this->infeasibility_measures[i]))) {
			dominated_entries++;
		}
	}
	return (dominated_entries <= this->number_dominated_entries); // point acceptable (dominated by <= this->number_dominated_entries)
}

/* check NON-MONOTONE sufficient reduction condition */
double NonmonotoneFilter::compute_actual_reduction(double current_objective, double current_residual, double trial_objective) {
	/* max penalty among most recent entries */
	double max_objective = current_objective;
	for (int i = 0; i < this->number_dominated_entries; i++) {
		double gamma;
		if (current_residual < this->infeasibility_measures[this->size-i]) {
			gamma = 1./this->constants.Gamma;
		}
		else {
			gamma = this->constants.Gamma;
		}
		double dash_objective = this->optimality_measures[this->size-i] + (this->infeasibility_measures[this->size-i] - current_residual)*gamma;
		max_objective = std::max(max_objective, dash_objective);
	}
	/* non-monotone actual reduction */
	return max_objective - trial_objective;
}
