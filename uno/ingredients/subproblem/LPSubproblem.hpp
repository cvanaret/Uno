// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSUBPROBLEM_H
#define UNO_LPSUBPROBLEM_H

#include "ActiveSetSubproblem.hpp"
#include "solvers/QP/LPSolver.hpp"
#include "tools/Options.hpp"

class LPSubproblem : public ActiveSetSubproblem {
public:
   LPSubproblem(size_t max_number_variables, size_t max_number_constraints, const Options& options);

   [[nodiscard]] Direction solve(Statistics& statistics, const ReformulatedProblem& problem, Iterate& current_iterate) override;
   [[nodiscard]] Direction compute_second_order_correction(const ReformulatedProblem& model, Iterate& trial_iterate) override;
   [[nodiscard]] PredictedOptimalityReductionModel generate_predicted_optimality_reduction_model(const ReformulatedProblem& problem, const Direction& direction) const override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] double get_proximal_coefficient() const override;

private:
   // pointer to allow polymorphism
   const std::unique_ptr<LPSolver> solver; /*!< Solver that solves the subproblem */

   void evaluate_functions(const ReformulatedProblem& problem, Iterate& current_iterate);
   [[nodiscard]] Direction solve_LP(const ReformulatedProblem& problem, Iterate& iterate);
};

#endif // UNO_LPSUBPROBLEM_H
