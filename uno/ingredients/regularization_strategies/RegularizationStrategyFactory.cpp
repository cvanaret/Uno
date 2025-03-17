// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "RegularizationStrategyFactory.hpp"
#include "PrimalRegularization.hpp"
#include "PrimalDualRegularization.hpp"
#include "NoRegularization.hpp"

namespace uno {
   std::unique_ptr<RegularizationStrategy<double>> RegularizationStrategyFactory::create(const std::string& strategy_name, const Options& options) {
      if (strategy_name == "primal") {
         return std::make_unique<PrimalRegularization<double>>(options);
      }
      else if (strategy_name == "primal_dual") {
         return std::make_unique<PrimalDualRegularization<double>>(options);
      }
      else if (strategy_name == "none") {
         return std::make_unique<NoRegularization<double>>();
      }
      throw std::invalid_argument("Convexification strategy " + strategy_name + " does not exist");
   }
} // namespace