// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSUBPROBLEM_H
#define UNO_QPSUBPROBLEM_H

#include "ActiveSetSubproblem.hpp"
#include "ingredients/subproblem/HessianModel.hpp"
#include "solvers/QP/QPSolver.hpp"
#include "tools/Options.hpp"

class QPSubproblem : public ActiveSetSubproblem {
public:
   QPSubproblem(size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros, const Options& options);

   [[nodiscard]] Direction solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate) override;
   [[nodiscard]] Direction compute_second_order_correction(const NonlinearProblem& model, Iterate& trial_iterate) override;
   [[nodiscard]] PredictedOptimalityReductionModel generate_predicted_optimality_reduction_model(const NonlinearProblem& problem,
         const Direction& direction) const override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] double get_proximal_coefficient() const override;

protected:
   // pointers to allow polymorphism
   const std::unique_ptr<HessianModel> hessian_model; /*!< Strategy to evaluate or approximate the Hessian */
   const std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */
   const double proximal_coefficient;

   void evaluate_functions(const NonlinearProblem& problem, Iterate& current_iterate);
   [[nodiscard]] Direction solve_QP(const NonlinearProblem& problem, Iterate& iterate);
};

#endif // UNO_QPSUBPROBLEM_H
