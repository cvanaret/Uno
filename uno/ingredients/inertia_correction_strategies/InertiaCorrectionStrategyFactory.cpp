// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <string>
#include "InertiaCorrectionStrategyFactory.hpp"
#include "PrimalInertiaCorrection.hpp"
#include "PrimalDualInertiaCorrection.hpp"
#include "NoInertiaCorrection.hpp"
#include "options/Options.hpp"

namespace uno {
   std::unique_ptr<InertiaCorrectionStrategy<double>> InertiaCorrectionStrategyFactory::create(const Options& options) {
      const std::string& strategy_name = options.get_string("inertia_correction_strategy");
      if (strategy_name == "primal") {
         return std::make_unique<PrimalInertiaCorrection<double>>(options);
      }
      else if (strategy_name == "primal_dual") {
         return std::make_unique<PrimalDualInertiaCorrection<double>>(options);
      }
      else if (strategy_name == "none") {
         return std::make_unique<NoInertiaCorrection<double>>();
      }
      throw std::invalid_argument("Inertia correction strategy " + strategy_name + " does not exist");
   }
} // namespace