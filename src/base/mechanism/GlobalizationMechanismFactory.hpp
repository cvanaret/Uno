#ifndef GLOBALIZATIONMECHANISMFACTORY_H
#define GLOBALIZATIONMECHANISMFACTORY_H

#include <iostream>
#include <memory>
#include <map>
#include "GlobalizationMechanism.hpp"
#include "GlobalizationStrategy.hpp"
#include "Options.hpp"

class GlobalizationMechanismFactory {
public:
    static std::unique_ptr<GlobalizationMechanism> create(const std::string& mechanism_type, ConstraintRelaxationStrategy&
    constraint_relaxation_strategy, const Options& options);
};

#endif // GLOBALIZATIONMECHANISMFACTORY_H
