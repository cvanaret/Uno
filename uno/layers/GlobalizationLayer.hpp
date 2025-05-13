// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONLAYER_H
#define UNO_GLOBALIZATIONLAYER_H

#include <memory>
#include "ingredients/globalization_mechanisms/GlobalizationMechanism.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"

namespace uno {
   // forward declarations
   class Iterate;
   class Options;
   class Statistics;

   class GlobalizationLayer {
   public:
      std::unique_ptr<GlobalizationStrategy> strategy;
      std::unique_ptr<GlobalizationMechanism> mechanism;

      GlobalizationLayer(size_t number_constraints, const Options& options):
         strategy(GlobalizationStrategyFactory::create(number_constraints, options)),
         mechanism(GlobalizationMechanismFactory::create(options)) {
      }

      void initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) const {
         this->strategy->initialize(statistics, initial_iterate, options);
         this->mechanism->initialize(statistics, options);
      }
   };
} // namespace

#endif //UNO_GLOBALIZATIONLAYER_H