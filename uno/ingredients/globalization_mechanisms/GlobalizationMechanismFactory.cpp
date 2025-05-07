// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "GlobalizationMechanism.hpp"
#include "GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_mechanisms/TrustRegionStrategy.hpp"
#include "ingredients/globalization_mechanisms/BacktrackingLineSearch.hpp"
#include "options/Options.hpp"

namespace uno {
   std::unique_ptr<GlobalizationMechanism> GlobalizationMechanismFactory::create(size_t number_constraints, size_t number_bounds_constraints,
         const Options& options) {
      const std::string& mechanism_type = options.get_string("globalization_mechanism");
       if (mechanism_type == "TR") {
           return std::make_unique<TrustRegionStrategy>(number_constraints, number_bounds_constraints, options);
       }
       else if (mechanism_type == "LS") {
           return std::make_unique<BacktrackingLineSearch>(number_constraints, number_bounds_constraints, options);
       }
       throw std::invalid_argument("GlobalizationMechanism " + mechanism_type + " is not supported");
   }

   std::vector<std::string> GlobalizationMechanismFactory::available_strategies() {
      return {"TR", "LS"};
   }
} // namespace