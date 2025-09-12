// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGYFACTORY_H

#include <array>
#include <memory>

namespace uno {
   // forward declarations
   class ConstraintRelaxationStrategy;
   class Options;

   class ConstraintRelaxationStrategyFactory {
   public:
      static std::unique_ptr<ConstraintRelaxationStrategy> create(bool unconstrained_model, const Options& options);

      constexpr static std::array available_strategies{
         "feasibility_restoration"
      };
   };
} // namespace

#endif // UNO_CONSTRAINTRELAXATIONSTRATEGYFACTORY_H