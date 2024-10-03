// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONMECHANISMFACTORY_H
#define UNO_GLOBALIZATIONMECHANISMFACTORY_H

#include <memory>
#include <vector>
//#include "GlobalizationMechanism.hpp"
//#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
//#include "tools/Options.hpp"

namespace uno {
   // forward declarations
   class ConstraintRelaxationStrategy;
   class GlobalizationMechanism;
   class Options;

   class GlobalizationMechanismFactory {
   public:
      static std::unique_ptr<GlobalizationMechanism> create(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);
      static std::vector<std::string> available_strategies();
   };
} // namespace

#endif // UNO_GLOBALIZATIONMECHANISMFACTORY_H
