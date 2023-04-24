// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_mechanism/TrustRegionStrategy.hpp"
#include "ingredients/globalization_mechanism/BacktrackingLineSearch.hpp"

std::unique_ptr<GlobalizationMechanism> GlobalizationMechanismFactory::create(Statistics& statistics,
      ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options) {
   const std::string& mechanism_type = options.get_string("globalization_mechanism");
    if (mechanism_type == "TR") {
        return std::make_unique<TrustRegionStrategy>(statistics, constraint_relaxation_strategy, options);
    }
    else if (mechanism_type == "LS") {
        return std::make_unique<BacktrackingLineSearch>(statistics, constraint_relaxation_strategy, options);
    }
    throw std::invalid_argument("GlobalizationMechanism " + mechanism_type + " is not supported");
}

std::vector<std::string> GlobalizationMechanismFactory::available_strategies() {
   return {"TR", "LS"};
}