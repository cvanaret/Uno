// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONSTRATEGYFACTORY_H
#define UNO_GLOBALIZATIONSTRATEGYFACTORY_H

#include <memory>
#include <vector>

namespace uno {
   // forward declarations
   class GlobalizationStrategy;
   class Options;

   class GlobalizationStrategyFactory {
   public:
      static std::unique_ptr<GlobalizationStrategy> create(const std::string& strategy_type, const Options& options);
      static std::vector<std::string> available_strategies();
   };
} // namespace

#endif // UNO_GLOBALIZATIONSTRATEGYFACTORY_H
