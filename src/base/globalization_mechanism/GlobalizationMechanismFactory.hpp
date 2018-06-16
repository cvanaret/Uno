#ifndef GLOBALIZATIONMECHANISMFACTORY_H
#define GLOBALIZATIONMECHANISMFACTORY_H

#include <iostream>
#include <memory>
#include <map>
#include "GlobalizationMechanism.hpp"
#include "GlobalizationStrategy.hpp"

class GlobalizationMechanismFactory {
	public:
		static std::shared_ptr<GlobalizationMechanism> create(const std::string& type, GlobalizationStrategy& globalization_strategy,
			std::map<std::string, std::string> default_values);
};

#endif // GLOBALIZATIONMECHANISMFACTORY_H
