#ifndef SUBPROBLEMFACTORY_H
#define SUBPROBLEMFACTORY_H

#include <iostream>
#include <memory>
#include "Subproblem.hpp"
#include "QPSolver.hpp"
#include "HessianEvaluation.hpp"

class SubproblemFactory {
	public:
		static std::shared_ptr<Subproblem> create(const std::string& type, QPSolver& solver, HessianEvaluation& hessian_evaluation, std::map<std::string, std::string> default_values);
};

#endif // SUBPROBLEMFACTORY_H
