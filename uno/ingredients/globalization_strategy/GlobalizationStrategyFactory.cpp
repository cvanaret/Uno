// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "GlobalizationStrategyFactory.hpp"
#include "l1MeritFunction.hpp"
#include "filter_method/FletcherFilterMethod.hpp"
#include "filter_method/WaechterFilterMethod.hpp"

namespace uno {
   std::unique_ptr <GlobalizationStrategy> GlobalizationStrategyFactory::create(const std::string& strategy_type, const Options& options) {
      if (strategy_type == "l1_merit") {
         return std::make_unique<l1MeritFunction>(options);
      }
      else if (strategy_type == "fletcher_filter_method") {
         return std::make_unique<FletcherFilterMethod>(options);
      }
      else if (strategy_type == "waechter_filter_method") {
         return std::make_unique<WaechterFilterMethod>(options);
      }
      throw std::invalid_argument("GlobalizationStrategy " + strategy_type + " is not supported");
   }

   std::vector<std::string> GlobalizationStrategyFactory::available_strategies() {
      return {"l1_merit", "fletcher_filter_strategy", "waechter_filter_strategy"};
   }
} // namespace
