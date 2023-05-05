// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSUBPROBLEM_H
#define UNO_QPSUBPROBLEM_H

#include "ActiveSetSubproblem.hpp"
#include "ingredients/subproblem/HessianModel.hpp"
#include "solvers/QP/QPSolver.hpp"
#include "tools/Options.hpp"

class QPSubproblem : public ActiveSetSubproblem {
public:
   QPSubproblem(Statistics& statistics, size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros,
         const Options& options);

   void generate_initial_iterate(const NonlinearProblem& problem, Iterate& initial_iterate) override;
   [[nodiscard]] Direction solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
         const WarmstartInformation& warmstart_information) override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;

protected:
   const bool use_regularization;
   // pointers to allow polymorphism
   const std::unique_ptr<HessianModel> hessian_model; /*!< Strategy to evaluate or approximate the Hessian */
   const std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */

   void evaluate_functions(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
         const WarmstartInformation& warmstart_information);
   [[nodiscard]] Direction solve_QP(const NonlinearProblem& problem, Iterate& iterate, const WarmstartInformation& warmstart_information);
};

#endif // UNO_QPSUBPROBLEM_H
