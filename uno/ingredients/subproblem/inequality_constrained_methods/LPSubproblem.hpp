// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSUBPROBLEM_H
#define UNO_LPSUBPROBLEM_H

#include <memory>
#include "InequalityConstrainedMethod.hpp"
#include "linear_algebra/COOSymmetricMatrix.hpp"
#include "solvers/LP/LPSolver.hpp"

class LPSubproblem : public InequalityConstrainedMethod {
public:
   LPSubproblem(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
         const Options& options);

   void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
   void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
         Direction& direction, const WarmstartInformation& warmstart_information) override;
   [[nodiscard]] const SymmetricMatrix<double>& get_lagrangian_hessian() const override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;

private:
   // pointer to allow polymorphism
   const std::unique_ptr<LPSolver> solver; /*!< Solver that solves the subproblem */
   const COOSymmetricMatrix<double> zero_hessian;

   void evaluate_functions(const OptimizationProblem& problem, Iterate& current_iterate, const WarmstartInformation& warmstart_information);
};

#endif // UNO_LPSUBPROBLEM_H
