// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSUBPROBLEM_H
#define UNO_LPSUBPROBLEM_H

#include <memory>
#include "InequalityConstrainedMethod.hpp"

namespace uno {
   // forward reference
   class LPSolver;

   class LPSubproblem : public InequalityConstrainedMethod {
   public:
      explicit LPSubproblem(const Options& options);
      ~LPSubproblem() override;

      void initialize(const OptimizationProblem& first_reformulation) override;
      void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
            Direction& direction, HessianModel& hessian_model, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& primal_direction) const override;

   private:
      // pointer to allow polymorphism
      const std::unique_ptr<LPSolver> solver;
   };
} // namespace

#endif // UNO_LPSUBPROBLEM_H