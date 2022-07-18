// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_GLOBALIZATIONMECHANISMFACTORY_H
#define UNO_GLOBALIZATIONMECHANISMFACTORY_H

#include "GlobalizationMechanism.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategy.hpp"
#include "tools/Options.hpp"

class GlobalizationMechanismFactory {
public:
    static std::unique_ptr<GlobalizationMechanism> create(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);
};

#endif // UNO_GLOBALIZATIONMECHANISMFACTORY_H
