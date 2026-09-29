// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYHANDLINGMETHODFACTORY_H
#define UNO_INEQUALITYHANDLINGMETHODFACTORY_H

#include <memory>
#include <vector>
#include "options/Options.hpp"

namespace uno {
   // forward declarations
   class OptimizationProblem;
   class InequalityHandlingMethod;

   class InequalityHandlingMethodFactory {
      public:
         static std::unique_ptr<InequalityHandlingMethod> create(const OptimizationProblem& problem, bool uses_trust_region,
            double objective_multiplier, const Options& options, std::vector<OptionOverride>& option_overrides);

         static std::vector<std::string> available_strategies();
   };
} // namespace

#endif // UNO_INEQUALITYHANDLINGMETHODFACTORY_H