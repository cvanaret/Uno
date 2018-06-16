#ifndef GLOBALIZATIONSTRATEGYFACTORY_H
#define GLOBALIZATIONSTRATEGYFACTORY_H

#include <iostream>
#include <memory>
#include "GlobalizationStrategy.hpp"
#include "LocalApproximation.hpp"

class GlobalizationStrategyFactory {
	public:
		static std::shared_ptr<GlobalizationStrategy> create(const std::string& type, LocalApproximation& local_approximation,
			std::map<std::string, std::string> default_values);
};

#endif // GLOBALIZATIONSTRATEGYFACTORY_H
