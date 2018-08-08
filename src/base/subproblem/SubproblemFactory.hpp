#ifndef SUBPROBLEMFACTORY_H
#define SUBPROBLEMFACTORY_H

#include <iostream>
#include <memory>
#include "Subproblem.hpp"
#include "QPSolver.hpp"

class SubproblemFactory {
	public:
		static std::shared_ptr<Subproblem> create(const std::string& type, QPSolver& solver, std::map<std::string, std::string> default_values);
};

#endif // SUBPROBLEMFACTORY_H
