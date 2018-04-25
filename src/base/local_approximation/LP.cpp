#include "LP.hpp"

LP::LP(int number_variables, int number_constraints):
		variable_lb(number_variables), variable_ub(number_variables),
		constraint_lb(number_constraints), constraint_ub(number_constraints) {
	this->number_variables = number_variables;
	this->number_constraints = number_constraints;
}

std::ostream& operator<< (std::ostream &stream, LP& lp) {
	stream << "min g^T d\n";
	stream << " d\n";
	
	/* variables */
	stream << lp.number_variables << " variables\n";
	for (int i = 0; i < lp.number_variables; i ++) {
		stream << lp.variable_lb[i] << " <= " << "d" << i << " <= " << lp.variable_ub[i] << "\n";
	}
	
	/* objective */
	stream << "g = ";
	for (std::map<int,double>::iterator it = lp.objective.begin(); it != lp.objective.end(); it++) {
		int index = it->first;
		double value = it->second;
		if (it != lp.objective.begin()) {
			stream << " + ";
		}
		stream << value << "*d" << index;
	}
	stream << "\n";
	
	/* constraints */
	stream << lp.number_constraints << " constraints\n";
	for (int j = 0; j < lp.number_constraints; j++) {
		stream << lp.constraint_lb[j] << " <= ";
		for (std::map<int,double>::iterator it = lp.constraints[j].begin(); it != lp.constraints[j].end(); it++) {
			int index = it->first;
			double value = it->second;
			if (it != lp.constraints[j].begin()) {
				stream << " + ";
			}
			stream << value << "*d" << index;
		}
		stream << " <= " << lp.constraint_ub[j] << "\n";
	}
	
	return stream;
}
