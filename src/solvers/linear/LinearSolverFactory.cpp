#include "Vector.hpp"
#include "LinearSolverFactory.hpp"
#ifdef HAS_MA57
#include "MA57Solver.hpp"
#endif

std::unique_ptr<LinearSolver> LinearSolverFactory::create(const std::string& linear_solver_name) {
    std::vector<std::string> possible_solvers;
#ifdef HAS_MA57
    if (linear_solver_name == "MA57") {
        return std::make_unique<MA57Solver>();
    }
    possible_solvers.push_back("MA57");
#endif
//#ifdef HAS_PARDISO
//    if (linear_solver_name == "PARDISO") {
//        return std::make_unique<PardisoSolver>();
//    }
//    possible_solvers.push_back("PARDISO");
//#endif
    throw std::invalid_argument("LinearSolver name " + linear_solver_name + " does not exist.");
       // " The possible options are: " + join(possible_solvers, ", "));
}
