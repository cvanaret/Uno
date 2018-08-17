#include <iostream>
#include <ctime>
#include "Argonot.hpp"
#include "Iterate.hpp"
#include "Logger.hpp"

Argonot::Argonot(GlobalizationMechanism& globalization_mechanism, int max_iterations):
		globalization_mechanism(globalization_mechanism), max_iterations(max_iterations) {
}

Result Argonot::solve(Problem& problem, std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers) {
	std::clock_t c_start = std::clock();
	
	int major_iterations = 0, minor_iterations = 0;
	
	INFO << "Problem " << problem.name << "\n";
	INFO << problem.number_variables << " variables, " << problem.number_constraints << " constraints\n";

	/* evaluate the initial point */
	Iterate current_iterate(problem, x, bound_multipliers, constraint_multipliers);
	INFO << "Initial iterate\n" << current_iterate << "\n";
	
	/* use the evaluation of the current point to initialize the strategies */
	this->globalization_mechanism.initialize(problem, current_iterate);
	
	try {
		/* check for convergence */
		while (!this->termination_criterion(current_iterate.status, major_iterations)) {
			major_iterations++;
			DEBUG << "\n\t\tARGONOT iteration " << major_iterations << "\n";
			INFO << "major: " << major_iterations << "\t";

			/* update the current point */
			
			current_iterate = this->globalization_mechanism.compute_iterate(problem, current_iterate);
			minor_iterations += this->globalization_mechanism.number_iterations;
			
			INFO << "constraints: " << current_iterate.residual << "\tobjective: " << current_iterate.objective << "\t";
			INFO << "status: " << current_iterate.status << "\n";
			DEBUG << "Next iterate\n" << current_iterate;
		}
	}
	catch (std::out_of_range& exception) {
		ERROR << exception.what();
	}
	catch (std::invalid_argument& exception) {
		ERROR << exception.what();
	}
	std::clock_t c_end = std::clock();
	double cpu_time = (c_end-c_start) / (double) CLOCKS_PER_SEC;

	Result result = {current_iterate,
					major_iterations,
					cpu_time,
					problem.number_eval_objective,
					problem.number_eval_constraints,
					problem.number_eval_jacobian,
					problem.number_eval_hessian,
					this->globalization_mechanism.globalization_strategy.subproblem.number_subproblems_solved};
	return result;
}

bool Argonot::termination_criterion(OptimalityStatus current_status, int iteration) {
	return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

void Result::display() {
	std::cout << "\n";
	std::cout << "ARGONOT v1: optimization summary\n";
	std::cout << "==============================\n";
	
	std::cout << "Status:\t\t\t";
	if (this->solution.status == KKT_POINT) {
		std::cout << "Feasible KKT point found\n";
	}
	else if (this->solution.status == FJ_POINT) {
		std::cout << "Infeasible FJ point found\n";
	}
	else if (this->solution.status == FEASIBLE_SMALL_STEP) {
		std::cout << "Feasible small step\n";
	}
	else if (this->solution.status == INFEASIBLE_SMALL_STEP) {
		std::cout << "Infeasible small step\n";
	}
	else { // NOT_OPTIMAL
		std::cout << "Irregular termination\n";
	}

	std::cout << "Objective value:\t" << this->solution.objective << "\n";
	std::cout << "Constraint residual:\t" << this->solution.residual << "\n";
	std::cout << "KKT residual:\t\t" << this->solution.KKTerror << "\n";
	std::cout << "Complementarity:\t" << this->solution.complementarity_error << "\n";

	std::cout << "Primal solution:\t";
	print_vector(std::cout, this->solution.x);
	
	std::cout << "Bound multipliers:\t\t";
	print_vector(std::cout, this->solution.bound_multipliers);
	std::cout << "Constraint multipliers:\t\t";
	print_vector(std::cout, this->solution.constraint_multipliers);
	
	std::cout << "CPU time:\t\t" << this->cpu_time << "s\n";
	std::cout << "Iterations:\t\t" << this->iteration << "\n";
	std::cout << "Objective evaluations:\t\t" << this->objective_evaluations << "\n";
	std::cout << "Constraints evaluations:\t\t" << this->constraint_evaluations << "\n";
	std::cout << "Jacobian evaluations:\t\t" << this->jacobian_evaluations << "\n";
	std::cout << "Hessian evaluations:\t\t" << this->hessian_evaluations << "\n";
	std::cout << "Number of subproblems solved:\t\t" << this->number_subproblems_solved << "\n";
	return;
}
