// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_GLOBALIZATIONSTRATEGYFACTORY_H
#define UNO_GLOBALIZATIONSTRATEGYFACTORY_H

#include <memory>
#include "GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

class GlobalizationStrategyFactory {
public:
   static std::unique_ptr<GlobalizationStrategy> create(const std::string& strategy_type, const Options& options);
};

#endif // UNO_GLOBALIZATIONSTRATEGYFACTORY_H
