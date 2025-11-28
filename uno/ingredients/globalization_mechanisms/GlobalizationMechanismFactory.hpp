// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONMECHANISMFACTORY_H
#define UNO_GLOBALIZATIONMECHANISMFACTORY_H

#include <array>
#include <memory>
#include "GlobalizationMechanism.hpp"

namespace uno {
   // forward declarations
   class Model;
   class Options;

   class GlobalizationMechanismFactory {
   public:
      static std::unique_ptr<GlobalizationMechanism> create(const Model& model, const Options& options);

      constexpr static std::array available_strategies{"TR", "LS"};
   };
} // namespace

#endif // UNO_GLOBALIZATIONMECHANISMFACTORY_H
