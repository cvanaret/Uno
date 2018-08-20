#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <fstream>
#include "AMPLModel.hpp"
#include "QPSolverFactory.hpp"
#include "SubproblemFactory.hpp"
#include "GlobalizationStrategyFactory.hpp"
#include "GlobalizationMechanismFactory.hpp"
#include "Argonot.hpp"
#include "Logger.hpp"
#include "InteriorPoint.hpp"

void test_argonot(std::string problem_name, std::map<std::string, std::string> options) {
    AMPLModel problem = AMPLModel(problem_name);

    /* create the local solver */
    std::shared_ptr<QPSolver> solver = QPSolverFactory::create("BQPD", problem, options);

    /* create the subproblem strategy */
    std::shared_ptr<Subproblem> subproblem = SubproblemFactory::create(options["subproblem"], *solver, options);

    /* create the globalization strategy */
    std::shared_ptr<GlobalizationStrategy> strategy = GlobalizationStrategyFactory::create(options["strategy"], *subproblem, options);

    /* create the globalization mechanism */
    std::shared_ptr<GlobalizationMechanism> mechanism = GlobalizationMechanismFactory::create(options["mechanism"], *strategy, options);

    int max_iterations = stoi(options["max_iterations"]);
    Argonot argonot = Argonot(*mechanism, max_iterations);

    /* initial primal and dual points */
    std::vector<double> x = problem.primal_initial_solution();
    std::vector<double> constraint_multipliers = problem.dual_initial_solution();
    std::vector<double> bound_multipliers(problem.number_variables);

    Result result = argonot.solve(problem, x, bound_multipliers, constraint_multipliers);
    result.display();
    return;
}

std::map<std::string, std::string> get_command_options(int argc, char* argv[], std::map<std::string, std::string>& options) {
    /* build the (argument, value) map */
    int i = 1;
    while (i < argc - 1) {
        std::string argument_i = std::string(argv[i]);

        if (argument_i[0] == '-' && i < argc - 1) {
            std::string value_i = std::string(argv[i + 1]);
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

void test_interior_point() {
    /* test of hs015 */
    AMPLModel problem = AMPLModel("../ampl_models/hs015");

    std::vector<double> x = problem.primal_initial_solution();
    std::vector<double> constraint_multipliers = problem.dual_initial_solution();
    std::vector<double> bound_multipliers(problem.number_variables);
    Iterate current_iterate(problem, x, bound_multipliers, constraint_multipliers);

    double radius = INFINITY;
    InteriorPoint ipm;
    ipm.compute_measures(problem, current_iterate);
    ipm.initialize(problem, current_iterate, problem.number_variables, problem.number_constraints, radius < INFINITY);
    ipm.compute_optimality_step(problem, current_iterate, radius);

    return;
}

int main(int argc, char* argv[]) {
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
            std::string problem_name = std::string(argv[argc - 1]);
            test_argonot(problem_name, options);
        }
    }
    else {
        test_interior_point();
    }
    return 0;
}
