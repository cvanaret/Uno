// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDMETHOD_H
#define UNO_INEQUALITYCONSTRAINEDMETHOD_H

#include <memory>
#include "../InequalityHandlingMethod.hpp"

namespace uno {
   class InequalityConstrainedMethod : public InequalityHandlingMethod {
   public:
      InequalityConstrainedMethod() = default;
      ~InequalityConstrainedMethod() override = default;

      void check_problem(const OptimizationProblem& problem, bool uses_trust_region) override;
      void initialize_statistics(Statistics& statistics) override;
      [[nodiscard]] std::unique_ptr<OptimizationProblem> reformulate(const OptimizationProblem& problem,
         Parameterization& parameterization) override;
      [[nodiscard]] bool update_parameterization(Statistics& statistics, const OptimizationProblem& problem,
         const Iterate& current_iterate, Parameterization& parameterization) override;

      void initialize_feasibility_problem(Iterate& current_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate, Evaluations& evaluations) override;
      [[nodiscard]] double proximal_coefficient() const override;

      [[nodiscard]] std::string get_name() const override;
   };
} // namespace

#endif // UNO_INEQUALITYCONSTRAINEDMETHOD_H