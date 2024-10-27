// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iomanip>
#include "Result.hpp"
#include "TerminationStatus.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   void Result::print(bool print_primal_dual_solution) const {
      DISCRETE << "Status:\t\t\t\t\t";
      if (this->solution.status == TerminationStatus::FEASIBLE_KKT_POINT) {
         DISCRETE << "Converged with feasible KKT point\n";
      }
      else if (this->solution.status == TerminationStatus::FEASIBLE_FJ_POINT) {
         DISCRETE << "Converged with feasible FJ point\n";
      }
      else if (this->solution.status == TerminationStatus::INFEASIBLE_STATIONARY_POINT) {
         DISCRETE << "Converged with infeasible stationary point\n";
      }
      else if (this->solution.status == TerminationStatus::FEASIBLE_SMALL_STEP) {
         DISCRETE << "Terminated with feasible small step\n";
      }
      else if (this->solution.status == TerminationStatus::INFEASIBLE_SMALL_STEP) {
         DISCRETE << "Failed with infeasible small step\n";
      }
      else if (this->solution.status == TerminationStatus::UNBOUNDED) {
         DISCRETE << "Terminated with unbounded problem\n";
      }
      else {
         DISCRETE << "Failed with suboptimal point\n";
      }

      DISCRETE << "Objective value:\t\t\t" << std::defaultfloat << std::setprecision(7) << this->solution.evaluations.objective << '\n';
      DISCRETE << "Primal feasibility:\t\t\t" << this->solution.primal_feasibility << '\n';

      DISCRETE << "┌ Stationarity residual:\t\t" << this->solution.residuals.stationarity << '\n';
      DISCRETE << "└ Complementarity residual:\t\t" << this->solution.residuals.complementarity << '\n';

      DISCRETE << "┌ Feasibility stationarity residual:\t" << this->solution.residuals.stationarity << '\n';
      DISCRETE << "└ Feasibility complementarity residual:\t" << this->solution.residuals.complementarity << '\n';

      DISCRETE << "┌ Infeasibility measure:\t\t" << this->solution.progress.infeasibility << '\n';
      DISCRETE << "│ Objective measure:\t\t\t" << this->solution.progress.objective(1.) << '\n';
      DISCRETE << "└ Auxiliary measure:\t\t\t" << this->solution.progress.auxiliary << '\n';

      if (print_primal_dual_solution) {
         DISCRETE << "Primal solution:\t\t\t"; print_vector(DISCRETE, view(this->solution.primals, 0, this->number_variables));
         DISCRETE << "┌ Constraint multipliers:\t\t"; print_vector(DISCRETE, this->solution.multipliers.constraints);
         DISCRETE << "│ Lower bound multipliers:\t\t"; print_vector(DISCRETE, view(this->solution.multipliers.lower_bounds, 0,
               this->number_variables));
         DISCRETE << "└ Upper bound multipliers:\t\t"; print_vector(DISCRETE, view(this->solution.multipliers.upper_bounds, 0,
               this->number_variables));
         DISCRETE << "┌ Constraint feasibility multipliers:\t"; print_vector(DISCRETE, this->solution.feasibility_multipliers.constraints);
         DISCRETE << "│ Lower bound feasibility multipliers:\t"; print_vector(DISCRETE, view(this->solution.feasibility_multipliers.lower_bounds, 0,
               this->number_variables));
         DISCRETE << "└ Upper bound feasibility multipliers:\t"; print_vector(DISCRETE, view(this->solution.feasibility_multipliers.upper_bounds, 0,
               this->number_variables));
         DISCRETE << "Objective multiplier:\t\t\t" << this->solution.objective_multiplier << '\n';
      }

      DISCRETE << "CPU time:\t\t\t\t" << this->cpu_time << "s\n";
      DISCRETE << "Iterations:\t\t\t\t" << this->iteration << '\n';
      DISCRETE << "Objective evaluations:\t\t\t" << this->objective_evaluations << '\n';
      DISCRETE << "Constraints evaluations:\t\t" << this->constraint_evaluations << '\n';
      DISCRETE << "Objective gradient evaluations:\t\t" << this->objective_gradient_evaluations << '\n';
      DISCRETE << "Jacobian evaluations:\t\t\t" << this->jacobian_evaluations << '\n';
      DISCRETE << "Hessian evaluations:\t\t\t" << this->hessian_evaluations << '\n';
      DISCRETE << "Number of subproblems solved:\t\t" << this->number_subproblems_solved << '\n';
   }
} // namespace
