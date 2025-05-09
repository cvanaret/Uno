// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_REGULARIZATIONSTRATEGYFACTORY_H
#define UNO_REGULARIZATIONSTRATEGYFACTORY_H

#include <memory>
#include "RegularizationStrategy.hpp"

namespace uno {
   // forward declaration
   class Options;

   class RegularizationStrategyFactory {
   public:
      static std::unique_ptr<RegularizationStrategy<double>> create(const Options& options);
   };
} // namespace

#endif // UNO_REGULARIZATIONSTRATEGYFACTORY_H