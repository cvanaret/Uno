// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONSTRATEGYFACTORY_H
#define UNO_GLOBALIZATIONSTRATEGYFACTORY_H

#include <array>
#include <memory>

namespace uno {
   // forward declarations
   class GlobalizationStrategy;
   class Model;
   class Options;

   class GlobalizationStrategyFactory {
   public:
      static std::unique_ptr<GlobalizationStrategy> create(const Model& model, const Options& options);

      constexpr static std::array available_strategies{
         "l1_merit", "fletcher_filter_method", "waechter_filter_method", "funnel_method"
      };
   };
} // namespace

#endif // UNO_GLOBALIZATIONSTRATEGYFACTORY_H