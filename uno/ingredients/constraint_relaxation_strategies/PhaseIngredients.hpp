// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PHASEINGREDIENTS_H
#define UNO_PHASEINGREDIENTS_H

#include <memory>
#include <string>

namespace uno {
   // forward declarations
   class GlobalizationStrategy;
   class InequalityHandlingMethod;

   struct PhaseIngredientsView {
      const std::string& problem_name;
      const std::string& phase_label;
      const OptimizationProblem& problem;
      InequalityHandlingMethod& inequality_handling_method;
      GlobalizationStrategy& globalization_strategy;
   };

   template <typename Problem>
   class PhaseIngredients {
   public:
      std::string problem_name; // optimality or feasibility
      std::string phase_label; // OPT or FEAS
      Problem problem;
      std::unique_ptr<InequalityHandlingMethod> inequality_handling_method{};
      std::unique_ptr<GlobalizationStrategy> globalization_strategy{};

      template <typename... Args>
      explicit PhaseIngredients(std::string problem_name, std::string phase_label, Args&&... args):
         problem_name(std::move(problem_name)), phase_label(std::move(phase_label)), problem(std::forward<Args>(args)...) {
      }

      // create a view to erase the type of Problem
      [[nodiscard]] PhaseIngredientsView view() const {
         return {this->problem_name, this->phase_label, this->problem, *this->inequality_handling_method,
            *this->globalization_strategy};
      }
   };
} // namespace

#endif //UNO_PHASEINGREDIENTS_H