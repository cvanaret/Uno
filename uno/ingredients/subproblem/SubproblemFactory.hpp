// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEMFACTORY_H
#define UNO_SUBPROBLEMFACTORY_H

#include <memory>
#include "Subproblem.hpp"
#include "ingredients/constraint_relaxation_strategy/NonlinearProblem.hpp"
#include "tools/Options.hpp"

class SubproblemFactory {
	public:
		static std::unique_ptr<Subproblem> create(size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros,
            const Options& options);
};

#endif // UNO_SUBPROBLEMFACTORY_H
