// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Result.hpp"

void Result::print(bool print_primal_dual_solution) const {
   std::cout << "Status:\t\t\t\t\t";
   if (this->status == FEASIBLE_KKT_POINT) {
      std::cout << "Converged with feasible KKT point\n";
   }
   else if (this->status == FJ_POINT) {
      std::cout << "Converged with feasible FJ point\n";
   }
   else if (this->status == INFEASIBLE_KKT_POINT) {
      std::cout << "Converged with infeasible KKT point\n";
   }
   else if (this->status == FEASIBLE_SMALL_STEP) {
      std::cout << "Terminated with feasible small step\n";
   }
   else if (this->status == INFEASIBLE_SMALL_STEP) {
      std::cout << "Terminated with infeasible small step\n";
   }
   else {
      std::cout << "Failed with error\n";
   }

   std::cout << "Objective value:\t\t\t" << this->solution.evaluations.objective << '\n';

   std::cout << "Constraint violation:\t\t\t" << this->solution.residuals.infeasibility << '\n';
   std::cout << "Optimality stationarity error:\t\t" << this->solution.residuals.optimality_stationarity << '\n';
   std::cout << "Feasibility stationarity error:\t\t" << this->solution.residuals.feasibility_stationarity << '\n';
   std::cout << "Optimality complementarity error:\t" << this->solution.residuals.optimality_complementarity << '\n';
   std::cout << "Feasibility complementarity error:\t" << this->solution.residuals.feasibility_complementarity << '\n';

   std::cout << "Infeasibility measure:\t\t\t" << this->solution.progress.infeasibility << '\n';
   std::cout << "Scaled optimality measure:\t\t" << this->solution.progress.scaled_optimality(1.) << '\n';
   std::cout << "Unscaled optimality measure:\t\t" << this->solution.progress.unscaled_optimality << '\n';

   if (print_primal_dual_solution) {
      std::cout << "Primal solution:\t\t\t"; print_vector(std::cout, this->solution.primals);
      if (not this->solution.multipliers.constraints.empty()) {
         std::cout << "Constraint multipliers:\t\t\t"; print_vector(std::cout, this->solution.multipliers.constraints);
      }
      std::cout << "Lower bound multipliers:\t\t"; print_vector(std::cout, this->solution.multipliers.lower_bounds);
      std::cout << "Upper bound multipliers:\t\t"; print_vector(std::cout, this->solution.multipliers.upper_bounds);
      std::cout << "Objective multiplier:\t\t\t" << this->solution.multipliers.objective << '\n';
   }

   std::cout << "CPU time:\t\t\t\t" << this->cpu_time << "s\n";
   std::cout << "Iterations:\t\t\t\t"	 << this->iteration << '\n';
   std::cout << "Objective evaluations:\t\t\t" << this->objective_evaluations << '\n';
   std::cout << "Constraints evaluations:\t\t" << this->constraint_evaluations << '\n';
   std::cout << "Jacobian evaluations:\t\t\t" << this->jacobian_evaluations << '\n';
   std::cout << "Hessian evaluations:\t\t\t" << this->hessian_evaluations << '\n';
   std::cout << "Number of subproblems solved:\t\t" << this->number_subproblems_solved << '\n';
}
