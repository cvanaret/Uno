#ifndef SUBPROBLEMFACTORY_H
#define SUBPROBLEMFACTORY_H

#include <iostream>
#include <memory>
#include "Subproblem.hpp"

class SubproblemFactory {
	public:
		static std::unique_ptr<Subproblem> create(const Problem& problem, size_t number_variables, const std::string& subproblem_type, const
		std::map<std::string, std::string>& options, bool use_trust_region);
};

#endif // SUBPROBLEMFACTORY_H
