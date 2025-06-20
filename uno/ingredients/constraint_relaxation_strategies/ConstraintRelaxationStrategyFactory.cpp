// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "l1Relaxation.hpp"
#include "UnconstrainedStrategy.hpp"
#include "options/Options.hpp"

namespace uno {
   std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(size_t number_constraints,
         const Options& options) {
      // set unconstrained strategy automatically
      if (number_constraints == 0) {
         return std::make_unique<UnconstrainedStrategy>(options);
      }
      const std::string constraint_relaxation_type = options.get_string("constraint_relaxation_strategy");
      if (constraint_relaxation_type == "feasibility_restoration") {
         return std::make_unique<FeasibilityRestoration>(options);
      }
      else if (constraint_relaxation_type == "l1_relaxation") {
         return std::make_unique<l1Relaxation>(options);
      }
      throw std::invalid_argument("ConstraintRelaxationStrategy " + constraint_relaxation_type + " is not supported");
   }

   std::vector<std::string> ConstraintRelaxationStrategyFactory::available_strategies() {
      return {"feasibility_restoration", "l1_relaxation"};
   }
} // namespace