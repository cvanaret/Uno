#include <iostream>
#include <vector>
#include <map>
#include "AMPLModel.hpp"
#include "BQPDSolver.hpp"
#include "QPApproximation.hpp"
#include "FilterStrategy.hpp"
#include "PenaltyStrategy.hpp"
#include "TrustRegion.hpp"
#include "LineSearch.hpp"
#include "TrustLineSearch.hpp"
#include "Argonot.hpp"
#include "Logger.hpp"

void test_argonot(std::string file_name) {
	AMPLModel problem = AMPLModel(file_name);
	
	BQPDSolver solver = BQPDSolver(problem.hessian_column_start, problem.hessian_row_number);
	QPApproximation local_approximation = QPApproximation(solver);
	
	double tolerance = 1e-6;
	
	
	FilterConstants filter_constants = {0.999, 0.001};
	Filter filter_restoration(filter_constants);
	Filter filter_optimality(filter_constants);
	LocalSolutionConstants step_constants = {0.1, 0.999};
	Tolerances filter_tolerances = {1e2, 1.25};
	FilterStrategy globalization_strategy = FilterStrategy(local_approximation, filter_restoration, filter_optimality, step_constants, filter_tolerances, tolerance);
	
	//PenaltyStrategy globalization_strategy = PenaltyStrategy(local_approximation, tolerance);
	
	
	double radius = 10.;
	TrustRegion globalization_mechanism = TrustRegion(globalization_strategy, radius);
	
	//TrustLineSearch globalization_mechanism = TrustLineSearch(globalization_strategy, radius);
	//LineSearch globalization_mechanism = LineSearch(globalization_strategy);
	
	int max_iterations = 1001;
	Argonot argonot = Argonot(globalization_mechanism, max_iterations);
	
	/* initial primal and dual points */
	std::vector<double> x = problem.primal_initial_solution();
	std::vector<double> multipliers = problem.dual_initial_solution();
	
	Result result = argonot.solve(problem, x, multipliers);
	result.display();
	return;
}

void handle_options(int argc, char *argv[]) {
	std::map<std::string, std::string> options;
	
	/* build the (argument, value) map */
	int i = 2;
	while (i < argc) {
		std::string argument_i = std::string(argv[i]);
		
		if (argument_i[0] == '-' && i < argc-1) {
			std::string value_i = std::string(argv[i+1]);
			options[argument_i] = value_i;
			i += 2;
		}
		else {
			std::cout << "Argument " << argument_i << " was ignored\n"; 
			i++;
		}
	}
	
	/* display the options */
	try {
		//std::string accept = options.at("-accept");
		//std::string globalization = options.at("-globalization");
		//std::string direction = options.at("-direction");
		//std::string accept = options.at("-accept");
		
		//std::cout << "Globalization: " << globalization << "\n";
		//std::cout << "Direction: " << direction << "\n";
		//std::cout << "Accept: " << accept << "\n";
	}
	catch (std::out_of_range) {
		std::cout << "Error in the options\n";
	}
}

Level Logger::logger_level = INFO;

int main (int argc, char* argv[]) {
	handle_options(argc, argv);
	
	if (1 < argc) {
		if (std::string(argv[1]).compare("-v") == 0) {
			std::cout << "Welcome in Argonot\n";
		}
		else {
			Logger::logger_level = DEBUG;
			test_argonot(std::string(argv[1]));
		}
	}
	return 0;
}
