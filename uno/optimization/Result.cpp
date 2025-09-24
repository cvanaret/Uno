// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iomanip>
#include "Result.hpp"
#include "SolutionStatus.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"

namespace uno {
   void Result::print(bool print_primal_dual_solution) const {
      DISCRETE << "Optimization status:\t\t\t" << optimization_status_to_message(this->optimization_status) << '\n';
      DISCRETE << "Solution status:\t\t\t" << solution_status_to_message(this->solution_status) << '\n';

      DISCRETE << "Objective value:\t\t\t" << std::defaultfloat << std::setprecision(7) << this->solution_objective << '\n';
      DISCRETE << "Primal feasibility:\t\t\t" << this->solution_primal_feasibility << '\n';

      DISCRETE << "┌ Stationarity residual:\t\t" << this->solution_dual_feasibility << '\n';
      DISCRETE << "│ Primal feasibility:\t\t\t" << this->solution_primal_feasibility << '\n';
      DISCRETE << "└ Complementarity residual:\t\t" << this->solution_complementarity << '\n';

      if (print_primal_dual_solution) {
         DISCRETE << "Primal solution:\t\t\t"; print_vector(DISCRETE, view(this->primal_solution, 0, this->number_variables));
         DISCRETE << "┌ Constraint multipliers:\t\t"; print_vector(DISCRETE, this->constraint_dual_solution);
         DISCRETE << "│ Lower bound multipliers:\t\t"; print_vector(DISCRETE, view(this->lower_bound_dual_solution, 0,
               this->number_variables));
         DISCRETE << "└ Upper bound multipliers:\t\t"; print_vector(DISCRETE, view(this->upper_bound_dual_solution, 0,
               this->number_variables));
      }

      DISCRETE << "CPU time:\t\t\t\t" << this->cpu_time << "s\n";
      DISCRETE << "Iterations:\t\t\t\t" << this->number_iterations << '\n';
      DISCRETE << "Objective evaluations:\t\t\t" << this->number_objective_evaluations << '\n';
      DISCRETE << "Constraints evaluations:\t\t" << this->number_constraint_evaluations << '\n';
      DISCRETE << "Objective gradient evaluations:\t\t" << this->number_objective_gradient_evaluations << '\n';
      DISCRETE << "Jacobian evaluations:\t\t\t" << this->number_jacobian_evaluations << '\n';
      DISCRETE << "Hessian evaluations:\t\t\t" << this->number_hessian_evaluations << '\n';
      DISCRETE << "Number of subproblems solved:\t\t" << this->number_subproblems_solved << '\n';
   }
} // namespace