// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategyFactory.hpp"
#include "FeasibilityRestoration.hpp"
#include "l1Relaxation.hpp"

std::unique_ptr<ConstraintRelaxationStrategy> ConstraintRelaxationStrategyFactory::create(Statistics& statistics, const Model& model,
      const Options& options) {
   const std::string constraint_relaxation_type = options.get_string("constraint_relaxation_strategy");
   if (constraint_relaxation_type == "feasibility_restoration") {
      return std::make_unique<FeasibilityRestoration>(statistics, model, options);
   }
   else if (constraint_relaxation_type == "l1_relaxation") {
      return std::make_unique<l1Relaxation>(statistics, model, options);
   }
   throw std::invalid_argument("ConstraintRelaxationStrategy " + constraint_relaxation_type + " is not supported");
}

std::vector<std::string> ConstraintRelaxationStrategyFactory::available_strategies() {
   return {"feasibility_restoration", "l1_relaxation"};
}