// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include <stdexcept>
#include "GlobalizationStrategy.hpp"
#include "GlobalizationStrategyFactory.hpp"
#include "l1MeritFunction.hpp"
#include "switching_methods/filter_methods/FletcherFilterMethod.hpp"
#include "switching_methods/filter_methods/WaechterFilterMethod.hpp"
#include "switching_methods/funnel_methods/FunnelMethod.hpp"

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
      else if (strategy_type == "funnel_method") {
         return std::make_unique<FunnelMethod>(options);
      }
      throw std::invalid_argument("GlobalizationStrategy " + strategy_type + " is not supported");
   }

   std::vector<std::string> GlobalizationStrategyFactory::available_strategies() {
      return {"l1_merit", "fletcher_filter_method", "waechter_filter_method", "funnel_method"};
   }
}
