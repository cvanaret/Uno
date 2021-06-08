#ifndef GLOBALIZATIONSTRATEGYFACTORY_H
#define GLOBALIZATIONSTRATEGYFACTORY_H

#include <iostream>
#include <memory>
#include "GlobalizationStrategy.hpp"
#include "Subproblem.hpp"

class GlobalizationStrategyFactory {
public:
   static std::unique_ptr<GlobalizationStrategy> create(const std::string& type, FeasibilityStrategy& feasibility_strategy,
         Subproblem& subproblem, const std::map<std::string, std::string>& options);
};

#endif // GLOBALIZATIONSTRATEGYFACTORY_H
