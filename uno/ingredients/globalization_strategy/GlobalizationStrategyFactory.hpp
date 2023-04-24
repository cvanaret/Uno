// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONSTRATEGYFACTORY_H
#define UNO_GLOBALIZATIONSTRATEGYFACTORY_H

#include <memory>
#include "GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

class GlobalizationStrategyFactory {
public:
   static std::unique_ptr<GlobalizationStrategy> create(Statistics& statistics, const std::string& strategy_type, const Options& options);
   static std::vector<std::string> available_strategies();
};

#endif // UNO_GLOBALIZATIONSTRATEGYFACTORY_H
