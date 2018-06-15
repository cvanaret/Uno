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
	this->constraints.reserve(this->max_size);
	this->objective.reserve(this->max_size);
	return;
}

/*  add (cons,objf) to the filter */
void Filter::add(double cons, double objf) {
	/* remove any redundant entries from the filter */

	/* find position in filter WITHOUT margin */
	int start = 0;
	while (start < this->size && cons >= this->constraints[start]) {
		start++;
	}

	/* find redundant entries starting from start */
	int end = start;
	while (end < this->size && objf <= this->objective[end]) {
		end++;
	}
	--end;

	int new_size = this->size;
	/* remove entries [start:end] from filter */
	int length = end - start + 1;
	if (0 < length) {
		this->shift_left(start, length);
		new_size -= length;
	}

	/* check sufficient space available for new entry (remove last entry, if not) */
	if (this->max_size <= this->size) {
		/* set upper bound */
		this->upper_bound = this->constants.Beta*std::max(this->upper_bound, this->constraints[this->size-1]);
		/* create space in filter: remove entry this->size-1 */
		new_size--;
	}

	/* add new entry to the filter */

	/*  find position of cons in filter (new entry will be index) */
	int index = 0;
	while (index < new_size && cons >= this->constants.Beta*this->constraints[index]) {
		index++;
	}

	/* shift entries by one to right to make room for new entry */
	if (index < new_size) {
		int length = 1;
		this->size = new_size; // filter size changes above
		this->shift_right(index, length);
	}

	/* add new entry to filter in position index */
	this->constraints[index] = cons;
	this->objective[index] = objf;
	this->size = new_size + 1;

	return;
}

/* query: return true if (cons,objf) acceptable, false otherwise */
bool Filter::query(double cons, double objf) {
	/* check upper bound first */
	if (this->constants.Beta*this->upper_bound <= cons) {
		return false;
	}

	/* find position of cons in filter (new entry will be index) */
	int index = 0;
	while (index < this->size && cons >= this->constants.Beta*this->constraints[index]) {
		index++;
	}

	/* now check acceptability */
	if (index == 0)
		return true; // acceptable as left-most entry
	else if (objf <= this->objective[index-1] - this->constants.Gamma*cons) 
		return true; // point acceptable
	else 
		return false; // point rejected by entry [index-1]
}

//! query_current_iterate: check acceptable wrt current point 
bool Filter::query_current_iterate(double curc, double curf, double newc, double newf) {
	return !((newf > curf - this->constants.Gamma*newc) && (newc > this->constants.Beta*curc));
}

//! shift_left: shift filter to left from start to end
void Filter::shift_left(int start, int length) {
	for (int i = start; i < this->size-length; i++) {
		this->constraints[i] = this->constraints[i+length];
		this->objective[i] = this->objective[i+length];
	}
	return;
}

//! shift_right: shift filter to right from end to start
void Filter::shift_right(int start, int length) {
	for (int i = this->size; i > start; i--) {
		this->constraints[i] = this->constraints[i-length];
		this->objective[i] = this->objective[i-length];
	}
	return;
}

double Filter::compute_actual_reduction(double current_objective, double current_residual, double new_objective) {
	return current_objective - new_objective;
}

//! print: print the complete filter
std::ostream& operator<< (std::ostream &stream, Filter& filter) {
	stream << "    Current filter (constraint residual, objective):\n";
	for (int i = 0; i < filter.size; i++) {
		//fprintf(stdout, "\t%16.8g %16.8g\n", this->constraints[i], this->objective[i]);
		stream << "\t" << filter.constraints[i] << "\t" << filter.objective[i] << "\n";
	}
	//fprintf(stdout, "    Upper bound = %16.8g\n", this->upper_bound);
	stream << "    Upper bound = " << filter.upper_bound << "\n";
	return stream;
}

/* NonmonotoneFilter class*/

NonmonotoneFilter::NonmonotoneFilter(FilterConstants& constants, int number_dominated_entries):
		Filter(constants), number_dominated_entries(number_dominated_entries) {
}

//!  add (cons,objf) to the filter
void NonmonotoneFilter::add(double cons, double objf) {
	int new_size = this->size;
	int dominated_entries;
	/* find entries in filter that are dominated by this->number_dominated_entries other entries */
	for (int i = 0; i < this->size; i++) {
		/* check whether ith entry dominated by (cons,objf) */
		if ((this->objective[i] > objf) && (this->constraints[i] > cons)) {
			dominated_entries = 1;
		}
		else {
			dominated_entries = 0;
		}
		
		/* find other filter entries that dominate ith entry */
		for (int j = 0; j < this->size; j++) {
			if ((this->objective[i] > this->objective[j]) && (this->constraints[i] > this->constraints[j]))  {
				dominated_entries++;
			}
		}
		if (dominated_entries > this->number_dominated_entries) {
			/* remove this entry */
			this->shift_left(i, 1);
			new_size--;
		}
	}

	/* check sufficient space available */
	if (new_size > this->max_size - 1) {
		/* create space in filter: remove entry 1 (oldest entry) */
		this->shift_left(1, 1);
		new_size--;
	}

	/* add new entry to filter in position new_size */
	this->constraints[new_size] = cons;
	this->objective[new_size] = objf;
	this->size = new_size + 1;

	return;
}

//! query: check if (cons,objf) acceptable; return 1 if yes, 0 if not
bool NonmonotoneFilter::query(double cons, double objf) {
	/* check upper bound first */
	if (cons >= this->constants.Beta*this->upper_bound) {
		return 0;
	}
	
	/* check acceptability by counting how many entries dominate */
	int dominated_entries = 0;
	for (int i = 0; i < this->size; i++) {
		if (((objf > this->objective[i] - this->constants.Gamma*cons) && (cons >= this->constants.Beta*this->constraints[i])) ||
				((objf >= this->objective[i] - this->constants.Gamma*cons) && (cons > this->constants.Beta*this->constraints[i]))) {
			dominated_entries++;
		}
	}
	
	return (dominated_entries <= this->number_dominated_entries); // point acceptable (dominated by <= this->number_dominated_entries)
}

//! query_current_iterate: check acceptable wrt current point 
bool NonmonotoneFilter::query_current_iterate(double curc, double curf, double newc, double newf) {
	int dominated_entries;
	
	/* check acceptability wrt current point (non-monotone) */
	if ((newf > curf - this->constants.Gamma*newc) && (newc > this->constants.Beta*curc)) 
		dominated_entries = 1;
	else
		dominated_entries = 0;
		
	for (int i = 0; i < this->size; i++) {
		if (((newf > this->objective[i] - this->constants.Gamma*newc) &&  (newc >= this->constants.Beta*this->constraints[i])) ||
				((newf >= this->objective[i] - this->constants.Gamma*newc) &&  (newc > this->constants.Beta*this->constraints[i]))) {
			dominated_entries++;
		}
	}
	
	return (dominated_entries <= this->number_dominated_entries); // point acceptable (dominated by <= this->number_dominated_entries)
}

/* check NON-MONOTONE sufficient reduction condition */
double NonmonotoneFilter::compute_actual_reduction(double current_objective, double current_residual, double new_objective) {
	/* max penalty among most recent entries */
	double max_objective = current_objective;
	for (int i = 0; i < this->number_dominated_entries; i++) {
		double gamma;
		if (current_residual < this->constraints[this->size-i]) {
			gamma = 1./this->constants.Gamma;
		}
		else {
			gamma = this->constants.Gamma;
		}
		double dash_objective = this->objective[this->size-i] + (this->constraints[this->size-i] - current_residual)*gamma;
		max_objective = std::max(max_objective, dash_objective);
	}
	/* non-monotone actual reduction */
	return max_objective - new_objective;
}
