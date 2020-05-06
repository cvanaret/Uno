#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include "AMPLModel.hpp"
#include "QPSolverFactory.hpp"
#include "SubproblemFactory.hpp"
#include "GlobalizationStrategyFactory.hpp"
#include "GlobalizationMechanismFactory.hpp"
#include "Argonot.hpp"
#include "Logger.hpp"
#include "Sl1QP.hpp"
#include "BQPDSolver.hpp"

void run_argonot(std::string problem_name, std::map<std::string, std::string> options) {
    // generate Hessians with a Fortran indexing (starting at 1) that is supported by solvers
    int fortran_indexing = 1;
    AMPLModel problem = AMPLModel(problem_name, fortran_indexing);
    
    /* create the subproblem strategy */
    std::shared_ptr<Subproblem> subproblem = SubproblemFactory::create(problem, options["subproblem"], options);

    /* create the globalization strategy */
    std::shared_ptr<GlobalizationStrategy> strategy = GlobalizationStrategyFactory::create(options["strategy"], *subproblem, options);

    /* create the globalization mechanism */
    std::shared_ptr<GlobalizationMechanism> mechanism = GlobalizationMechanismFactory::create(options["mechanism"], *strategy, options);

    int max_iterations = stoi(options["max_iterations"]);
    Argonot argonot = Argonot(*mechanism, max_iterations);

    /* initial primal and dual points */
    //std::cout << "Getting the initial point\n";
    std::vector<double> x = problem.primal_initial_solution();
    Multipliers multipliers(problem.number_variables, problem.number_constraints);
    multipliers.constraints = problem.dual_initial_solution();

    Result result = argonot.solve(problem, x, multipliers);
    /* remove auxiliary variables */
    result.solution.x.resize(problem.number_variables);
    result.solution.multipliers.lower_bounds.resize(problem.number_variables);
    result.solution.multipliers.upper_bounds.resize(problem.number_variables);
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
    } catch (std::out_of_range) {
    }
    return;
}

std::map<std::string, std::string> get_default_values(std::string file_name) {
    std::ifstream file;
    file.open(file_name);

    /* register the default values */
    std::map<std::string, std::string> default_values;
    std::string key, value;
    std::string line;
    while (std::getline(file, line)) {
        if (line != "" && line.find("#") != 0) {
            std::istringstream iss;
            iss.str(line);
            iss >> key >> value;
            std::cout << "Default value " << key << " = " << value << "\n";
            default_values[key] = value;
        }
    }
    file.close();
    return default_values;
}

//void test_sl1QP() {
//    AMPLModel problem = AMPLModel("../ampl_models/hs015");
//    BQPDSolver solver(problem.hessian_column_start, problem.hessian_row_number);
//    Sl1QP sl1qp(solver);
//
//    std::vector<double> x = problem.primal_initial_solution();
//    std::cout << "Original x is:     ";
//    print_vector(std::cout, x);
//    Multipliers multipliers(problem.number_variables, problem.number_constraints);
//    multipliers.constraints = problem.dual_initial_solution();
//
//    Iterate current_iterate = sl1qp.initialize(problem, x, multipliers, 2, false);
//    solver.allocate(current_iterate.x.size(), problem.number_constraints);
//
//    std::cout << "Reformulated x is: ";
//    print_vector(std::cout, current_iterate.x);
//
//    double trust_region_radius = INFINITY;
//    SubproblemSolution solution = sl1qp.compute_optimality_step(problem, current_iterate, trust_region_radius);
//    std::cout << solution << "\n";
//
//    double c1 = dot(solution.x, current_iterate.constraints_jacobian[0]);
//    double c2 = dot(solution.x, current_iterate.constraints_jacobian[1]);
//    std::cout << "Value of constraints: " << c1 << " and " << c2 << "\n";
//}

void test_matrix() {
    int fortran_indexing = 0;
    ArgonotMatrix matrix(3, fortran_indexing);
    matrix.add_term(1., 0, 0);
    matrix.add_term(2., 0, 2);
    matrix.add_term(3., 1, 2);
    matrix.add_term(4., 2, 0);
    matrix.add_term(5., 2, 1);
    matrix.add_term(6., 2, 2);
    CSCMatrix csc_matrix = matrix.to_CSC();
    std::cout << "Before\n" << csc_matrix;
    
    csc_matrix = csc_matrix.add_identity_multiple(100.);
    
    std::cout << "After\n" << csc_matrix;
}

void test_factorization() {
    int fortran_indexing = 1;
    int n = 2;
    COOMatrix coo_matrix(n, fortran_indexing);
    coo_matrix.add_term(6050.0001, 0, 0);
    coo_matrix.add_term(-2774, 0, 1);
    coo_matrix.add_term(1e-4, 1, 1);
    
    MA57Solver solver;
    MA57Factorization factorization = solver.factorize(coo_matrix);
    std::cout << "Dimension: " << coo_matrix.dimension << "\n";
    std::cout << "Singular ? " << factorization.matrix_is_singular() << "\n";
    std::cout << "Rank ? " << factorization.rank() << "\n";
    std::cout << "Negative eigenvalues ? " << factorization.number_negative_eigenvalues() << "\n";
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
            /* run Argonot */
            run_argonot(problem_name, options);
        }
    }
    return 0;
}
