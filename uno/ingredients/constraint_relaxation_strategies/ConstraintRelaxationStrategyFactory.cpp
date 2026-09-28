// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "NoRelaxation.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(const Model& model,
         bool use_trust_region, const Options& options, std::vector<OptionOverride>& option_overrides) {
      const std::string constraint_relaxation_type = options.get_string("constraint_relaxation_strategy");
      // figure out whether there are constraints altogether
      if (model.number_constraints == 0) {
         DEBUG << "The model is unconstrained, picking no relaxation\n";
         option_overrides.emplace_back("constraint_relaxation_strategy", constraint_relaxation_type, "no_relaxation",
            "unconstrained problem");
         return std::make_unique<NoRelaxation>(model, options);
      }
      // from now on, there are constraints
      if (constraint_relaxation_type == "feasibility_restoration") {
         DEBUG << "Picking feasibility restoration\n";
         return std::make_unique<FeasibilityRestoration>(model, use_trust_region, options, option_overrides);
      }
      throw std::invalid_argument("ConstraintRelaxationStrategy " + constraint_relaxation_type + " is not supported");
   }
} // namespace