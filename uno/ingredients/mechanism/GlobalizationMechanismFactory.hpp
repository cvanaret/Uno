#ifndef GLOBALIZATIONMECHANISMFACTORY_H
#define GLOBALIZATIONMECHANISMFACTORY_H

#include "GlobalizationMechanism.hpp"
#include "ingredients/constraint_relaxation/ConstraintRelaxationStrategy.hpp"
#include "tools/Options.hpp"

class GlobalizationMechanismFactory {
public:
    static std::unique_ptr<GlobalizationMechanism> create(std::string_view mechanism_type, ConstraintRelaxationStrategy&
    constraint_relaxation_strategy, const Options& options);
};

#endif // GLOBALIZATIONMECHANISMFACTORY_H
