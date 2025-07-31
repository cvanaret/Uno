// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "UnconstrainedStrategy.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(size_t number_constraints,
         const Options& options) {
      // set unconstrained strategy automatically
      if (number_constraints == 0) {
         INFO << "The model is unconstrained, picking an unconstrained constraint relaxation strategy\n";
         return std::make_unique<UnconstrainedStrategy>(options);
      }
      const std::string constraint_relaxation_type = options.get_string("constraint_relaxation_strategy");
      if (constraint_relaxation_type == "feasibility_restoration") {
         return std::make_unique<FeasibilityRestoration>(options);
      }
      throw std::invalid_argument("ConstraintRelaxationStrategy " + constraint_relaxation_type + " is not supported");
   }
} // namespace