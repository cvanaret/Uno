#ifndef SUBPROBLEMFACTORY_H
#define SUBPROBLEMFACTORY_H

#include <memory>
#include "Subproblem.hpp"
#include "tools/Options.hpp"

class SubproblemFactory {
	public:
		static std::unique_ptr<Subproblem> create(const Problem& problem, size_t max_number_variables, const Options& options);
};

#endif // SUBPROBLEMFACTORY_H
