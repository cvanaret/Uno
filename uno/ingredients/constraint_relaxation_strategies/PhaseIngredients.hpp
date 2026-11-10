// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PHASEINGREDIENTS_H
#define UNO_PHASEINGREDIENTS_H

#include <memory>

namespace uno {
   class GlobalizationStrategy;
   class InequalityHandlingMethod;

   struct PhaseIngredientsView {
      const OptimizationProblem& problem;
      InequalityHandlingMethod& inequality_handling_method;
      GlobalizationStrategy& globalization_strategy;
   };

   template <typename Problem>
   class PhaseIngredients {
   public:
      Problem problem;
      std::unique_ptr<InequalityHandlingMethod> inequality_handling_method{};
      std::unique_ptr<GlobalizationStrategy> globalization_strategy{};

      template <typename... Args>
      PhaseIngredients(Args&&... args): problem(std::forward<Args>(args)...) { }

      [[nodiscard]] PhaseIngredientsView view() const {
         return {this->problem, *this->inequality_handling_method, *this->globalization_strategy};
      }
   };
} // namespace

#endif //UNO_PHASEINGREDIENTS_H