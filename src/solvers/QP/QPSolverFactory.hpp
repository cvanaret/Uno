#ifndef QPSOLVERFACTORY_H
#define QPSOLVERFACTORY_H

#include <iostream>
#include <memory>
#include <map>
#include "QPSolver.hpp"
#include "Problem.hpp"

class QPSolverFactory {
	public:
		static std::shared_ptr<QPSolver> create(const std::string& name, Problem& problem, std::map<std::string, std::string> default_values);
};

#endif // QPSOLVERFACTORY_H
