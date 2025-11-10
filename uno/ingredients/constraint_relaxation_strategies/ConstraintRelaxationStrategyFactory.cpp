// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "NoRelaxation.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(const Model& model,
         bool use_trust_region, const Options& options) {
      // figure out whether we need to relax constraints altogether
      if (model.number_constraints == 0) {
         INFO << "The model is unconstrained, picking no relaxation\n";
         return std::make_unique<NoRelaxation>(model, options);
      }
      const std::string constraint_relaxation_type = options.get_string("constraint_relaxation_strategy");
      if (constraint_relaxation_type == "feasibility_restoration") {
         return std::make_unique<FeasibilityRestoration>(model, use_trust_region, options);
      }
      throw std::invalid_argument("ConstraintRelaxationStrategy " + constraint_relaxation_type + " is not supported");
   }
} // namespace