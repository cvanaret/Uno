#ifndef QPSOLVERFACTORY_H
#define QPSOLVERFACTORY_H

#include <iostream>
#include <memory>
#include "QPSolver.hpp"
#include "Problem.hpp"

class QPSolverFactory {
	public:
		static std::shared_ptr<QPSolver> create(const std::string& name, Problem& problem);
};

#endif // QPSOLVERFACTORY_H
