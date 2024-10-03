// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGYFACTORY_H

#include <memory>
#include <vector>

namespace uno {
   // forward declarations
   class ConstraintRelaxationStrategy;
   class Model;
   class Options;

   class ConstraintRelaxationStrategyFactory {
   public:
      static std::unique_ptr<ConstraintRelaxationStrategy> create(const Model& model, const Options& options);
      static std::vector<std::string> available_strategies();
   };
} // namespace

#endif // UNO_CONSTRAINTRELAXATIONSTRATEGYFACTORY_H
