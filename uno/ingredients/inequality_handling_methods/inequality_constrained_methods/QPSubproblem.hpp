// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSUBPROBLEM_H
#define UNO_QPSUBPROBLEM_H

#include <memory>
#include "InequalityConstrainedMethod.hpp"

namespace uno {
   // forward reference
   class QPSolver;

   class QPSubproblem : public InequalityConstrainedMethod {
   public:
      explicit QPSubproblem(const Options& options);
      ~QPSubproblem() override;

      void initialize(const OptimizationProblem& first_reformulation, const HessianModel& hessian_model) override;
      void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
         Direction& direction, HessianModel& hessian_model, double trust_region_radius, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& vector) const override;

      [[nodiscard]] std::string get_strategy_combination() const override;

   protected:
      const bool enforce_linear_constraints_at_initial_iterate;
      // pointer to allow polymorphism
      const std::unique_ptr<QPSolver> solver;
   };
} // namespace

#endif // UNO_QPSUBPROBLEM_H