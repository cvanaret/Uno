// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONMECHANISMFACTORY_H
#define UNO_GLOBALIZATIONMECHANISMFACTORY_H

#include "GlobalizationMechanism.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategy.hpp"
#include "tools/Options.hpp"

class GlobalizationMechanismFactory {
public:
   static std::unique_ptr<GlobalizationMechanism> create(Statistics& statistics, ConstraintRelaxationStrategy& constraint_relaxation_strategy,
         const Options& options);
   static std::vector<std::string> available_strategies();
};

#endif // UNO_GLOBALIZATIONMECHANISMFACTORY_H
