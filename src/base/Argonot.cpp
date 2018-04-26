#include <iostream>
#include <ctime>
#include "Argonot.hpp"
#include "Iterate.hpp"
#include "Logger.hpp"

Argonot::Argonot(GlobalizationMechanism& globalization_mechanism, int max_iterations):
		globalization_mechanism(globalization_mechanism) {
	this->max_iterations = max_iterations;
}

Result Argonot::solve(Problem& problem, std::vector<double>& x, std::vector<double>& multipliers) {
	std::clock_t c_start = std::clock();
	
	int major_iterations = 0;
	int minor_iterations = 0;
	
	INFO << "Problem " << problem.name << "\n";
	INFO << problem.number_variables << " variables, " << problem.number_constraints << " constraints\n";

	/* evaluate the initial point */
	Iterate current_point(problem, x, multipliers);
	INFO << "Initial point\n" << current_point << "\n";
	
	/* use the evaluation of the current point to initialize the strategies */
	this->globalization_mechanism.initialize(problem, current_point);
	
	/* check for convergence */
	while (!this->termination_criterion(current_point.status, major_iterations)) {
		major_iterations++;
		DEBUG << "\n\t\tARGONOT iteration " << major_iterations << "\n";

		/* update the current point */
		INFO << "major: " << major_iterations << "\t";
		current_point = this->globalization_mechanism.compute_iterate(problem, current_point);
		INFO << "constraints: " << current_point.residual << "\tobjective: " << current_point.objective << "\t";
		INFO << "status: " << current_point.status << "\n";
		DEBUG << "Next point\n" << current_point;
		
		minor_iterations += this->globalization_mechanism.number_iterations;
	}
	
	std::clock_t c_end = std::clock();
	double cpu_time = (c_end-c_start) / (double) CLOCKS_PER_SEC;

	Result result = {current_point,
					major_iterations,
					cpu_time,
					problem.number_eval_objective,
					problem.number_eval_constraints,
					problem.number_eval_hessian,
					this->globalization_mechanism.globalization_strategy.local_approximation.number_subproblems_solved};
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
		std::cout << "Argonot reached max iter\n";
	}

	std::cout << "Objective value:\t" << this->solution.objective << "\n";
	std::cout << "Constraint residual:\t" << this->solution.residual << "\n";
	std::cout << "KKT residual:\t\t" << this->solution.KKTerror << "\n";
	std::cout << "Complementarity:\t" << this->solution.complementarity_error << "\n";

	std::cout << "Primal solution:\t";
	print_vector(std::cout, this->solution.x);
	
	std::cout << "Dual solution:\t\t";
	print_vector(std::cout, this->solution.multipliers);
	
	std::cout << "CPU time:\t\t" << this->cpu_time << "s\n";
	std::cout << "Iterations:\t\t" << this->iteration << "\n";
	std::cout << "Objective evaluations:\t\t" << this->objective_evaluations << "\n";
	std::cout << "Constraints evaluations:\t\t" << this->constraint_evaluations << "\n";
	std::cout << "Hessian evaluations:\t\t" << this->hessian_evaluations << "\n";
	std::cout << "Number of subproblems solved:\t\t" << this->number_subproblems_solved << "\n";
	return;
}

/* check termination conditions */
	//ifail = CheckIfail(ifail, optimal, nlp->phase, MajorIter, radius, *h);

////! check ifail & return appropriate error code
//int CheckIfail(int ifail, int optimal, Phase phase, int MajorIter, double rho, double h) {
  //if (ifail == 0){
    //if      ((optimal == KKT_POINT) && (phase == OPTIMALITY)) return(0); // optimal solution
    //else if ((optimal == KKT_POINT) && (phase == RESTORATION)) return(3); // min. of infeasibility
    //else if (optimal == SMALL_STEP)                   return(2); // step got too small
    //else if (MajorIter >= IterLim.major)     return(6); // iteration limit
    //else if (rho <= tolerance)                return(5); // TR got too small 
    //else if ((h <= tolerance)&&(phase == RESTORATION))  return(4); // h<eps, but QP inf.
  //}
  //else if ((ifail == 1)&&(MajorIter == 0))   return(11);// IEEE error in obj/con
  //else if ((ifail == 2)&&(MajorIter == 0))   return(12);// IEEE error in grad/hess
  //else if (ifail == 1)                       return(7); // IEEE error in obj/con
  //else if (ifail == 2)                       return(7); // IEEE error in grad/hess
  //else if (ifail >= 4)                       return(8); // crash in QP solver
  //fprintf(stdout,"Unknown Error in CheckIfail:\n");
  //fprintf(stdout,"   ifail = %d ; optimal = %d ; phase = %d; iter = %d\n",
	  //ifail, optimal, phase, MajorIter);
  //fprintf(stdout,"   rho = %16.8g ; h = %16.8g\n", rho, h);
  //return(13); // unknown error
//}
