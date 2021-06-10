#ifndef GLOBALIZATIONMECHANISMFACTORY_H
#define GLOBALIZATIONMECHANISMFACTORY_H

#include <iostream>
#include <memory>
#include <map>
#include "GlobalizationMechanism.hpp"
#include "GlobalizationStrategy.hpp"

class GlobalizationMechanismFactory {
public:
    static std::unique_ptr<GlobalizationMechanism> create(const std::string& mechanism_type, ConstraintRelaxationStrategy&
    constraint_relaxation_strategy, const std::map<std::string, std::string>& options);
};

#endif // GLOBALIZATIONMECHANISMFACTORY_H
