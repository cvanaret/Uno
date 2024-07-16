// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSUBPROBLEM_H
#define UNO_QPSUBPROBLEM_H

#include "InequalityConstrainedMethod.hpp"
#include "ingredients/subproblem/HessianModel.hpp"
#include "solvers/QP/QPSolver.hpp"

class QPSubproblem : public InequalityConstrainedMethod {
public:
   QPSubproblem(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
         size_t number_hessian_nonzeros, const Options& options);

   void initialize_statistics(Statistics& statistics, const Options& options) override;
   void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
   void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
         Direction& direction, const WarmstartInformation& warmstart_information) override;
   [[nodiscard]] const SymmetricMatrix<double>& get_lagrangian_hessian() const override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;

protected:
   const bool use_regularization;
   const bool enforce_linear_constraints_at_initial_iterate;
   // pointers to allow polymorphism
   const std::unique_ptr<HessianModel> hessian_model; /*!< Strategy to evaluate or approximate the Hessian */
   const std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */

   void evaluate_functions(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
         const WarmstartInformation& warmstart_information);
};

#endif // UNO_QPSUBPROBLEM_H
