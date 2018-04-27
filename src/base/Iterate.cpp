#include "Iterate.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

Iterate::Iterate(Problem& problem, std::vector<double>& x, std::vector<double>& multipliers):
		x(x), multipliers(multipliers) {
	
	/* project point onto bounds */
	for(int i = 0; i < problem.number_variables; i++) {
		this->x[i] = std::max(this->x[i], problem.variable_lb[i]);
		this->x[i] = std::min(this->x[i], problem.variable_ub[i]);
	}

	/* objective */
	this->objective = problem.objective(this->x);

	/* constraints */
	this->constraints = problem.evaluate_constraints(this->x);
	this->residual = problem.l1_inf_norm(this->constraints);
	this->infeasibility_measure = this->residual;
	this->feasibility_measure = this->objective;
	
	/* Jacobian and Hessian will be evaluated only when necessary */
	this->is_objective_gradient_computed = false;
	this->is_constraint_jacobian_computed = false;
	this->is_hessian_computed = false;
	this->is_hessian_computed = false;
	
	/* status */
	this->status = NOT_OPTIMAL;
}

void Iterate::set_objective_gradient(std::map<int,double>& objective_gradient) {
	this->objective_gradient = objective_gradient;
	this->is_objective_gradient_computed = true;
	return;
}

void Iterate::compute_constraint_jacobian(Problem& problem) {
	if (!this->is_constraint_jacobian_computed) {
		std::vector<std::map<int,double> > constraint_jacobian(problem.number_constraints);
		for (int j = 0; j < problem.number_constraints; j++) {
			constraint_jacobian[j] = problem.constraint_sparse_gradient(j, this->x);
		}
		this->constraint_jacobian = constraint_jacobian;
		this->is_constraint_jacobian_computed = true;
	}
	return;
}

void Iterate::compute_hessian(Problem& problem, double objective_multiplier, std::vector<double>& constraint_multipliers) {
	if (!this->is_hessian_computed) {
		DEBUG << "\t\tMultipliers: " << objective_multiplier << " | ";
		print_vector(DEBUG, constraint_multipliers, 20);
		
		this->hessian = problem.lagrangian_hessian(this->x, objective_multiplier, constraint_multipliers);
		this->is_hessian_computed = true;
	}
	return;
}

std::ostream& operator<< (std::ostream &stream, Iterate& iterate) {
	stream << "Primal: ";
	print_vector(stream, iterate.x, 50);
	
	stream << "Dual: ";
	print_vector(stream, iterate.multipliers, 50);
	
	stream << "Objective: " << iterate.objective << "\n";
	
	//stream << "Constraints:";
	//for (double cj: iterate.constraints) {
	//	stream << " " << cj;
	//}
	//stream << "\n";
	
	stream << "Residual: " << iterate.residual << "\n";
	
	return stream;
}

std::ostream& operator<< (std::ostream &stream, OptimalityStatus& status) {
	if (status == NOT_OPTIMAL) {
		stream << "not optimal";
	}
	else if (status == KKT_POINT) {
		stream << "KKT point";
	}
	else if (status == FJ_POINT) {
		stream << "FJ point";
	}
	else if (status == FEASIBLE_SMALL_STEP) {
		stream << "feasible small step";
	}
	else if (status == INFEASIBLE_SMALL_STEP) {
		stream << "infeasible small step";
	}
	return stream;
}
