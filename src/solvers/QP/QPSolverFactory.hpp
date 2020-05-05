#ifndef QPSOLVERFACTORY_H
#define QPSOLVERFACTORY_H

#include <memory>
#include "QPSolver.hpp"
#include "Problem.hpp"

class QPSolverFactory {
	public:
		static std::shared_ptr<QPSolver> create(const std::string& QP_solver, int number_variables, int number_constraints, int maximum_number_nonzeros);
};

#endif // QPSOLVERFACTORY_H
