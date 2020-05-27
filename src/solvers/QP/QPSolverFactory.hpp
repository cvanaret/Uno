#ifndef QPSOLVERFACTORY_H
#define QPSOLVERFACTORY_H

#include <memory>
#include "QPSolver.hpp"
#include "Problem.hpp"

class QPSolverFactory {
	public:
		static std::shared_ptr<QPSolver> create(const std::string& QP_solver_name, int number_variables, int number_constraints, int maximum_number_nonzeros, bool quadratic_programming);
};

#endif // QPSOLVERFACTORY_H
