#include "Utils.hpp"
#include "LinearSolverFactory.hpp"
#ifdef MA57
#include "MA57Solver.hpp"
#endif

std::shared_ptr<LinearSolver> LinearSolverFactory::create(const std::string& linear_solver_name, int number_variables, int number_constraints, int maximum_number_nonzeros) {
    std::vector<std::string> possible_solvers;
#ifdef MA57
    if (linear_solver_name == "MA57") {
        return std::make_shared<MA57Solver>();
    }
    possible_solvers.push_back("MA57");
#endif
//#ifdef PARDISO
//    if (linear_solver_name == "PARDISO") {
//        return std::make_shared<PardisoSolver>();
//    }
//    possible_solvers.push_back("PARDISO");
//#endif
    throw std::invalid_argument("LinearSolver name " + linear_solver_name + " does not exist. The possible options are: " + join(possible_solvers, ", "));
}
