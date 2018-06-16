#include "QPSolverFactory.hpp"
#include "BQPDSolver.hpp"

std::shared_ptr<QPSolver> QPSolverFactory::create(const std::string& name, Problem& problem, std::map<std::string, std::string> default_values) {
	if (name == "BQPD") {
		return std::make_shared<BQPDSolver>(problem.hessian_column_start, problem.hessian_row_number);
	}
	else {
		throw std::invalid_argument("QPSolver name " + name + " does not exist");
	}
}
