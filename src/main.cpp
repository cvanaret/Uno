#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <fstream>
#include "AMPLModel.hpp"
#include "QPSolverFactory.hpp"
#include "LocalApproximationFactory.hpp"
#include "GlobalizationStrategyFactory.hpp"
#include "GlobalizationMechanismFactory.hpp"
#include "Argonot.hpp"
#include "Logger.hpp"
#include "MA57Solver.hpp"

void test_argonot(std::string problem_name, std::map<std::string, std::string> options) {
	AMPLModel problem = AMPLModel(problem_name);
	
	/* create the local solver */
	std::shared_ptr<QPSolver> solver = QPSolverFactory::create("BQPD", problem, options);
	
	/* create the local approximation strategy */
	std::shared_ptr<LocalApproximation> local_approximation = LocalApproximationFactory::create(options["local_approximation"], *solver, options);
	
	/* create the globalization strategy */
	std::shared_ptr<GlobalizationStrategy> strategy = GlobalizationStrategyFactory::create(options["strategy"], *local_approximation, options);
	
	/* create the globalization mechanism */
	std::shared_ptr<GlobalizationMechanism> mechanism = GlobalizationMechanismFactory::create(options["mechanism"], *strategy, options);
	
	int max_iterations = stoi(options["max_iterations"]);
	Argonot argonot = Argonot(*mechanism, max_iterations);
	
	/* initial primal and dual points */
	std::vector<double> x = problem.primal_initial_solution();
	std::vector<double> multipliers = problem.dual_initial_solution();

	Result result = argonot.solve(problem, x, multipliers);
	result.display();
	return;
}

std::map<std::string, std::string> get_command_options(int argc, char* argv[], std::map<std::string, std::string>& options) {
	/* build the (argument, value) map */
	int i = 1;
	while (i < argc-1) {
		std::string argument_i = std::string(argv[i]);

		if (argument_i[0] == '-' && i < argc-1) {
			std::string value_i = std::string(argv[i+1]);
			std::cout << "(" << argument_i << ", " << value_i << ")\n";
			options[argument_i.substr(1)] = value_i;
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
		std::string logger_level = options.at("logger");
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

std::map<std::string, std::string> get_default_values(std::string file_name) {
	std::ifstream file;
	file.open(file_name);
	
	std::map<std::string, std::string> default_values;
	std::string key, value;
	while (file >> key >> value) {
		std::cout << "Default value " << key << " = " << value << "\n";
		default_values[key] = value;
	}
	file.close();
	return default_values;
}

int main (int argc, char* argv[]) {
	if (1 < argc) {
		/* get the default values */
		std::map<std::string, std::string> options = get_default_values("argonot.cfg");
		/* get the command line options */
		get_command_options(argc, argv, options);
	
		set_logger(options);
		
		if (std::string(argv[1]).compare("-v") == 0) {
			std::cout << "Welcome in Argonot\n";
		}
		else {
			std::string problem_name = std::string(argv[argc-1]);
			test_argonot(problem_name, options);
		}
	}
	else {
		/* test of hs015 */
		std::vector<double> x = {0.49, 2.1};
		std::vector<double> lambda = {700., 0.};
		double mu = 1e-7;
		std::vector<double> s = {x[0]*x[1] - 1., x[0] + x[1]*x[1]};
		double v = 0.5 - x[0];
		
		int n = 8;
		COOMatrix matrix(n, 0);
		
		matrix.add_term(1200.*x[0]*x[0] + 2. - 400.*x[1], 0, 0);
		matrix.add_term(-400.*x[0] + lambda[0], 0, 1);
		matrix.add_term(200. + 2.*lambda[1], 1, 1);
		
		matrix.add_term(x[1], 0, 5);
		matrix.add_term(x[0], 1, 5);
		matrix.add_term(1., 0, 6);
		matrix.add_term(2.*x[1], 1, 6);
		
		matrix.add_term(1., 0, 7);
		
		matrix.add_term(mu/(s[0]*s[0]), 2, 2);
		matrix.add_term(mu/(s[1]*s[1]), 3, 3);
		
		matrix.add_term(1., 2, 5);
		matrix.add_term(1., 3, 6);
		
		matrix.add_term(mu/(v*v), 4, 4);
		
		matrix.add_term(-1., 4, 7);
		
		/* right hand side */
		std::vector<double> rhs = {
			-400.*x[0]*x[0]*x[0] - 2.*x[0] + 400.*x[0]*x[1] + 2.,
			-200.*x[1] + 200*x[0]*x[0],
			mu/s[0],
			mu/s[1],
			mu/v,
			-x[0]*x[1] + 1. + s[0],
			-x[0] - x[1]*x[1] + s[1],
			-x[0] + 0.5 - v
		};
		
		MA57Solver solver;
		LocalSolution solution = solver.solve(matrix, rhs);
		
		std::cout << "solution";
		for (unsigned int k = 0; k < solution.x.size(); k++) {
			std::cout << " " << solution.x[k];
		}
		std::cout << "\n";
	
	}
	return 0;
}
