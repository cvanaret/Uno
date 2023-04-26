// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEMFACTORY_H
#define UNO_SUBPROBLEMFACTORY_H

#include <memory>
#include "Subproblem.hpp"
#include "tools/Options.hpp"

class SubproblemFactory {
	public:
		static std::unique_ptr<Subproblem> create(Statistics& statistics, size_t max_number_variables, size_t max_number_constraints,
            size_t max_number_jacobian_nonzeros, size_t max_number_hessian_nonzeros, const Options& options);

      static std::vector<std::string> available_strategies();
};

#endif // UNO_SUBPROBLEMFACTORY_H
