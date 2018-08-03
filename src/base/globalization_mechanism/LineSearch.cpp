#include <cmath>
#include "LineSearch.hpp"
#include "Logger.hpp"

LineSearch::LineSearch(GlobalizationStrategy& globalization_strategy, int max_iterations, double ratio):
		GlobalizationMechanism(globalization_strategy, max_iterations), ratio(ratio) {
}

Iterate LineSearch::compute_iterate(Problem& problem, Iterate& current_iterate) {
	/* compute a trial direction */
	LocalSolution solution = this->globalization_strategy.compute_step(problem, current_iterate, INFINITY);
	
	/* step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ... */
	this->step_length = 1.;
	bool is_accepted = false;
	this->number_iterations = 0;
	
	while (!this->termination(is_accepted, this->number_iterations)) {
		this->number_iterations++;
		DEBUG << "\n\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << this->step_length << "\n";
			
		try {
			/* check whether the trial step is accepted */
			is_accepted = this->globalization_strategy.check_step(problem, current_iterate, solution, this->step_length);
		}
		catch (const std::invalid_argument& e) {
			WARNING << RED << e.what() << RESET "\n";
			is_accepted = false;
		}
		
		if (is_accepted) {
			/* print summary */
			DEBUG << CYAN "LS trial point accepted\n" RESET;
			INFO << "minor: " << std::fixed << this->number_iterations << "\t";
			INFO << "step length: " << step_length << "\t";
			// TODO: if stragegy == penalty, the step norm has no meaning!
			INFO << "step norm: " << std::fixed << (step_length*solution.norm) << "\t";
		}
		else {
			/* decrease the step length */
			this->step_length *= this->ratio;
		}
	}
	return current_iterate;
}

bool LineSearch::termination(bool is_accepted, int iteration) {
	if (is_accepted) {
		return true;
	}
	else if (this->max_iterations < iteration) {
		throw std::out_of_range("Line-search iteration limit reached");
	}
	/* step length gets too small */
	if (this->step_length < 1e-16) {
		throw std::out_of_range("Step length became too small");
	}
	return false;
}

/*
 * Interpolation functions
 */

double LineSearch::quadratic_interpolation(Problem& problem, Iterate& current_iterate, std::vector<double> direction, double steplength) {
	std::cout << "Current point: "; print_vector(std::cout, current_iterate.x);
	std::cout << "Direction: "; print_vector(std::cout, direction);
	
	/* compute trial point */
	std::vector<double> trial_point(problem.number_variables);
	for (int i = 0; i < problem.number_variables; i++) {
		trial_point[i] = current_iterate.x[i] + steplength*direction[i];
	}
	/* evaluate trial point */
	double phi_alpha0 = problem.objective(trial_point);
	std::cout << "phi(alpha0) = f(x + alpha0*p) = " << phi_alpha0 << "\n";
	/* compute dot product */
	if (!current_iterate.is_objective_gradient_computed) {
		std::map<int,double> objective_gradient = problem.objective_sparse_gradient(current_iterate.x);
		current_iterate.set_objective_gradient(objective_gradient);
	}
	double phi_prime_0 = dot(direction, current_iterate.objective_gradient);
	std::cout << "phi'(0) = nabla f(x)^T p = " << phi_prime_0 << "\n";
	
	/* compute the minimum of the quadratic */
	double a = (phi_alpha0 - current_iterate.objective - phi_prime_0*steplength)/(steplength*steplength);
	double b = phi_prime_0;
	std::cout << "a = " << a << ", b = " << b << "\n";
	return this->minimize_quadratic(a, b);
}

double LineSearch::cubic_interpolation(Problem& problem, Iterate& current_iterate, std::vector<double> direction, double steplength1, double steplength2) {
	std::cout << "Current point: "; print_vector(std::cout, current_iterate.x);
	std::cout << "Direction: "; print_vector(std::cout, direction);
	
	/* compute trial points */
	std::vector<double> trial_point1(problem.number_variables);
	std::vector<double> trial_point2(problem.number_variables);
	for (int i = 0; i < problem.number_variables; i++) {
		trial_point1[i] = current_iterate.x[i] + steplength1*direction[i];
		trial_point2[i] = current_iterate.x[i] + steplength2*direction[i];
	}
	/* evaluate trial points */
	double phi_alpha1 = problem.objective(trial_point1);
	double phi_alpha2 = problem.objective(trial_point2);
	std::cout << "phi(alpha1) = f(x + alpha1*p) = " << phi_alpha1 << "\n";
	std::cout << "phi(alpha2) = f(x + alpha2*p) = " << phi_alpha2 << "\n";
	/* compute dot product */
	if (!current_iterate.is_objective_gradient_computed) {
		std::map<int,double> objective_gradient = problem.objective_sparse_gradient(current_iterate.x);
		current_iterate.set_objective_gradient(objective_gradient);
	}
	double phi_prime_0 = dot(direction, current_iterate.objective_gradient);
	std::cout << "phi'(0) = nabla f(x)^T p = " << phi_prime_0 << "\n";
	
	/* compute the minimum of the cubic */
	double det = steplength1*steplength1*steplength2*steplength2*(steplength1 - steplength2);
	std::cout << "Det = " << det << "\n";
	double K1 = phi_alpha1 - current_iterate.objective - steplength1*phi_prime_0;
	double K2 = phi_alpha2 - current_iterate.objective - steplength2*phi_prime_0;
	
	double a = (steplength2*steplength2*K1 - steplength1*steplength1*K2)/det;
	double b = (-steplength2*steplength2*steplength2*K1 + steplength1*steplength1*steplength1*K2)/det;//phi_prime_0;
	double c = phi_prime_0;
	std::cout << "a = " << a << ", b = " << b << ", c = " << c << "\n";
	if (a == 0.) {
		return this->minimize_quadratic(b, c);
	}
	else {
		return this->minimize_cubic(a, b, c);
	}
}

/* return the minimum of x -> ax^2 + bx + R */
double LineSearch::minimize_quadratic(double a, double b) {
	return -b/(2.*a);
}

/* return the minimum of x -> ax^3 + bx^2 + cx + R */
double LineSearch::minimize_cubic(double a, double b, double c) {
	return (-b + std::sqrt(b*b - 3.*a*c))/(3.*a);
}
