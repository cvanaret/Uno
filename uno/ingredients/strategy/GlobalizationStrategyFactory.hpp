#ifndef GLOBALIZATIONSTRATEGYFACTORY_H
#define GLOBALIZATIONSTRATEGYFACTORY_H

#include <iostream>
#include <memory>
#include "GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

class GlobalizationStrategyFactory {
public:
   static std::unique_ptr<GlobalizationStrategy> create(const std::string& strategy_type, const Options& options);
};

#endif // GLOBALIZATIONSTRATEGYFACTORY_H
