// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INERTIACORRECTIONSTRATEGYFACTORY_H
#define UNO_INERTIACORRECTIONSTRATEGYFACTORY_H

#include <memory>
#include "InertiaCorrectionStrategy.hpp"

namespace uno {
   // forward declaration
   class Options;

   class InertiaCorrectionStrategyFactory {
   public:
      static std::unique_ptr<InertiaCorrectionStrategy<double>> create(const Options& options);
   };
} // namespace

#endif // UNO_INERTIACORRECTIONSTRATEGYFACTORY_H