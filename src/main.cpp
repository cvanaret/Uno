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

TrustRegion compute_tr(GlobalizationStrategy& globalization_strategy) {
	double radius = 10.;
	TrustRegion globalization_mechanism = TrustRegion(globalization_strategy, radius);
	return globalization_mechanism;
}

void test_argonot(std::string file_name, std::map<std::string, std::string> options) {
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

std::map<std::string, std::string> get_command_options(int argc, char *argv[]) {
	std::map<std::string, std::string> options;
	
	/* build the (argument, value) map */
	int i = 1;
	while (i < argc-1) {
		std::string argument_i = std::string(argv[i]);

		if (argument_i[0] == '-' && i < argc-1) {
			std::string value_i = std::string(argv[i+1]);
			std::cout << "(" << argument_i << ", " << value_i << ")\n";
			options[argument_i] = value_i;
			i += 2;
		}
		else {
			std::cout << "Argument " << argument_i << " was ignored\n"; 
			i++;
		}
	}
	
	return options;
}

Level Logger::logger_level = INFO;

void set_logger(std::map<std::string, std::string> options) {
	Logger::logger_level = INFO;
	
	try {
		std::string logger_level = options.at("-logger");
		if (logger_level.compare("ERROR") == 0) {
			Logger::logger_level = ERROR;
		}
		else if (logger_level.compare("WARNING") == 0) {
			Logger::logger_level = WARNING;
		}
		else if (logger_level.compare("INFO") == 0) {
			Logger::logger_level = INFO;
		}
		else if (logger_level.compare("DEBUG") == 0) {
			Logger::logger_level = DEBUG;
		}
		else if (logger_level.compare("DEBUG1") == 0) {
			Logger::logger_level = DEBUG1;
		}
		else if (logger_level.compare("DEBUG2") == 0) {
			Logger::logger_level = DEBUG2;
		}
		else if (logger_level.compare("DEBUG3") == 0) {
			Logger::logger_level = DEBUG3;
		}
		else if (logger_level.compare("DEBUG4") == 0) {
			Logger::logger_level = DEBUG4;
		}
	}
	catch (std::out_of_range) {
	}
	return;
}

int main (int argc, char* argv[]) {
	if (1 < argc) {
		/* get the command line options */
		std::map<std::string, std::string> options = get_command_options(argc, argv);
	
		set_logger(options);
		
		if (std::string(argv[1]).compare("-v") == 0) {
			std::cout << "Welcome in Argonot\n";
		}
		else {
			test_argonot(std::string(argv[argc-1]), options);
		}
	}
	return 0;
}
