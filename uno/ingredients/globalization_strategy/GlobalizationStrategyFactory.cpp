// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationStrategyFactory.hpp"
#include "MeritFunction.hpp"
#include "filter_strategy/LeyfferFilterStrategy.hpp"
#include "filter_strategy/WaechterFilterStrategy.hpp"

std::unique_ptr <GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& strategy_type, const Options& options) {
   if (strategy_type == "merit") {
      return std::make_unique<MeritFunction>(options);
   }
   else if (strategy_type == "leyffer-filter-strategy") {
      return std::make_unique<LeyfferFilterStrategy>(options);
   }
   else if (strategy_type == "waechter-filter-strategy") {
      return std::make_unique<WaechterFilterStrategy>(options);
   }
   throw std::invalid_argument("GlobalizationStrategy " + strategy_type + " is not supported");
}

std::vector<std::string> GlobalizationStrategyFactory::available_strategies() {
   return {"merit", "leyffer-filter-strategy", "waechter-filter-strategy"};
}