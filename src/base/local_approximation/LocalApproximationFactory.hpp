#ifndef LOCALAPPROXIMATIONFACTORY_H
#define LOCALAPPROXIMATIONFACTORY_H

#include <iostream>
#include <memory>
#include "LocalApproximation.hpp"
#include "QPSolver.hpp"

class LocalApproximationFactory {
	public:
		static std::shared_ptr<LocalApproximation> create(const std::string& type, QPSolver& solver, std::map<std::string, std::string> default_values);
};

#endif // LOCALAPPROXIMATIONFACTORY_H
