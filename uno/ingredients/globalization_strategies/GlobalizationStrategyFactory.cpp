// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include <stdexcept>
#include "GlobalizationStrategy.hpp"
#include "GlobalizationStrategyFactory.hpp"
#include "l1MeritFunction.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "switching_methods/filter_methods/FletcherFilterMethod.hpp"
#include "switching_methods/filter_methods/WaechterFilterMethod.hpp"
#include "switching_methods/funnel_methods/FunnelMethod.hpp"
#include "tools/Logger.hpp"

namespace uno {
   std::unique_ptr <GlobalizationStrategy> GlobalizationStrategyFactory::create(const Model& model, const Options& options) {
      // set unconstrained strategy automatically
      if (model.number_constraints == 0) {
         INFO << "The model is unconstrained, picking a merit function as globalization strategy\n";
         return std::make_unique<l1MeritFunction>(options);
      }
      const std::string& strategy_type = options.get_string("globalization_strategy");
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
}