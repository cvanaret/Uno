// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPEQPSUBPROBLEM_H
#define UNO_LPEQPSUBPROBLEM_H

#include "ActiveSetSubproblem.hpp"
#include "ingredients/subproblem/HessianModel.hpp"
#include "solvers/QP/QPSolver.hpp"
#include "tools/Options.hpp"

class LPEQPSubproblem : public ActiveSetSubproblem {
public:
   LPEQPSubproblem(size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros, const Options& options);

   void initialize(Statistics& statistics, const NonlinearProblem& problem, Iterate& first_iterate) override;
   [[nodiscard]] Direction solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate) override;
   [[nodiscard]] Direction compute_second_order_correction(const NonlinearProblem& model, Iterate& trial_iterate) override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;

protected:
   // pointers to allow polymorphism
   const std::unique_ptr<HessianModel> hessian_model; /*!< Strategy to evaluate or approximate the Hessian */
   const std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */

   void evaluate_functions(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate);
   [[nodiscard]] Direction solve_QP(const NonlinearProblem& problem, Iterate& iterate);
   [[nodiscard]] Direction solve_LP(const NonlinearProblem& problem, Iterate& iterate);
};

#endif // UNO_LPEQPSUBPROBLEM_H
