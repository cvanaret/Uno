// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSUBPROBLEM_H
#define UNO_LPSUBPROBLEM_H

#include <memory>
#include "InequalityConstrainedMethod.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"

namespace uno {
   // forward reference
   class LPSolver;

   class LPSubproblem : public InequalityConstrainedMethod {
   public:
      LPSubproblem(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
            const Options& options);
      ~LPSubproblem();

      void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
            Direction& direction, const WarmstartInformation& warmstart_information) override;

   private:
      // pointer to allow polymorphism
      const std::unique_ptr<LPSolver> solver; /*!< Solver that solves the subproblem */
      const SymmetricMatrix<size_t, double> zero_hessian;

      void evaluate_functions(const OptimizationProblem& problem, Iterate& current_iterate, const WarmstartInformation& warmstart_information);
   };
} // namespace

#endif // UNO_LPSUBPROBLEM_H