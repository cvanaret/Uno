#ifndef UNO_SUBPROBLEMFACTORY_H
#define UNO_SUBPROBLEMFACTORY_H

#include <memory>
#include "Subproblem.hpp"
#include "ingredients/constraint_relaxation/NonlinearProblem.hpp"
#include "tools/Options.hpp"

class SubproblemFactory {
	public:
		static std::unique_ptr<Subproblem> create(const NonlinearProblem& problem, const Options& options);
};

#endif // UNO_SUBPROBLEMFACTORY_H
