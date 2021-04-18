#ifndef SUBPROBLEMFACTORY_H
#define SUBPROBLEMFACTORY_H

#include <iostream>
#include <memory>
#include "Subproblem.hpp"

class SubproblemFactory {
	public:
		static std::unique_ptr<Subproblem> create(Problem& problem, const std::string& type, std::map<std::string, std::string> default_values, bool use_trust_region, bool scale_residuals);
};

#endif // SUBPROBLEMFACTORY_H
