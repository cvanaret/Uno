#ifndef SUBPROBLEMFACTORY_H
#define SUBPROBLEMFACTORY_H

#include <iostream>
#include <memory>
#include "Subproblem.hpp"
#include "Options.hpp"

class SubproblemFactory {
	public:
		static std::unique_ptr<Subproblem> create(const Problem& problem, size_t number_variables, const std::string& subproblem_type, const
		Options& options, bool use_trust_region);
};

#endif // SUBPROBLEMFACTORY_H
