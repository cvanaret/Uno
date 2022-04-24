#ifndef UNO_SUBPROBLEMFACTORY_H
#define UNO_SUBPROBLEMFACTORY_H

#include <memory>
#include "Subproblem.hpp"
#include "ingredients/constraint_relaxation/ReformulatedProblem.hpp"
#include "tools/Options.hpp"

class SubproblemFactory {
	public:
		static std::unique_ptr<Subproblem> create(const ReformulatedProblem& problem, const Options& options);
};

#endif // UNO_SUBPROBLEMFACTORY_H
