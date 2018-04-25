#include "QP.hpp"

QP::QP(int number_variables, int number_constraints, const Matrix& hessian):
		variable_lb(number_variables), variable_ub(number_variables),
		constraint_lb(number_constraints), constraint_ub(number_constraints),
		constraints(number_constraints),
		hessian(hessian) {
	this->number_variables = number_variables;
	this->number_constraints = number_constraints;
}

std::ostream& operator<< (std::ostream &stream, QP& qp) {
	stream << "min 1/2 d^T W d + g^T d\n";
	stream << " d\n";
	
	/* variables */
	stream << qp.number_variables << " variables\n";
	for (int i = 0; i < qp.number_variables; i ++) {
		stream << qp.variable_lb[i] << " <= " << "d" << i << " <= " << qp.variable_ub[i] << "\n";
	}
	
	/* objective */
	stream << "g = ";
	for (std::map<int,double>::iterator it = qp.objective.begin(); it != qp.objective.end(); it++) {
		int index = it->first;
		double value = it->second;
		if (it != qp.objective.begin()) {
			stream << " + ";
		}
		stream << value << "*d" << index;
	}
	stream << "\n";
	
	/* constraints */
	stream << qp.number_constraints << " constraints\n";
	for (int j = 0; j < qp.number_constraints; j++) {
		stream << qp.constraint_lb[j] << " <= ";
		for (std::map<int,double>::iterator it = qp.constraints[j].begin(); it != qp.constraints[j].end(); it++) {
			int index = it->first;
			double value = it->second;
			if (it != qp.constraints[j].begin()) {
				stream << " + ";
			}
			stream << value << "*d" << index;
		}
		stream << " <= " << qp.constraint_ub[j] << "\n";
	}
	
	/* Hessian */
	stream << qp.hessian;
	
	return stream;
}
