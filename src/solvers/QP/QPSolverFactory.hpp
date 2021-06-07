#ifndef QPSOLVERFACTORY_H
#define QPSOLVERFACTORY_H

#include <memory>
#include "QPSolver.hpp"
#include "Problem.hpp"

class QPSolverFactory {
	public:
		static std::unique_ptr<QPSolver> create(const std::string& QP_solver_name, size_t number_variables, size_t number_constraints,
		      size_t maximum_number_nonzeros, bool quadratic_programming);
};

#endif // QPSOLVERFACTORY_H
