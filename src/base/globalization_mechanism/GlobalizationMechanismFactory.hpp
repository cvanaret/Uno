#ifndef GLOBALIZATIONMECHANISMFACTORY_H
#define GLOBALIZATIONMECHANISMFACTORY_H

#include <iostream>
#include <memory>
#include "GlobalizationMechanism.hpp"
#include "GlobalizationStrategy.hpp"

class GlobalizationMechanismFactory {
	public:
		static std::shared_ptr<GlobalizationMechanism> create(const std::string& type, GlobalizationStrategy& globalization_strategy);
};

#endif // GLOBALIZATIONMECHANISMFACTORY_H
