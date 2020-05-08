#ifndef SUBPROBLEMFACTORY_H
#define SUBPROBLEMFACTORY_H

#include <iostream>
#include <memory>
#include "Subproblem.hpp"
#include "HessianEvaluation.hpp"

class SubproblemFactory {
	public:
		static std::shared_ptr<Subproblem> create(Problem& problem, const std::string& type, std::map<std::string, std::string> default_values, bool use_trust_region);
};

#endif // SUBPROBLEMFACTORY_H
