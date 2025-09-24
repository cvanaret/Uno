// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RESULT_H
#define UNO_RESULT_H

#include "Iterate.hpp"
#include "OptimizationStatus.hpp"

namespace uno {
   struct Result {
      Result() = delete;

      const size_t number_variables;
      const size_t number_constraints;
      const OptimizationStatus optimization_status;
      const SolutionStatus solution_status;
      const double solution_objective;
      const double solution_primal_feasibility;
      const double solution_dual_feasibility;
      const double solution_complementarity;
      Vector<double> primal_solution;
      Vector<double> constraint_dual_solution;
      Vector<double> lower_bound_dual_solution;
      Vector<double> upper_bound_dual_solution;
      const size_t number_iterations;
      const double cpu_time;
      const size_t number_objective_evaluations;
      const size_t number_constraint_evaluations;
      const size_t number_objective_gradient_evaluations;
      const size_t number_jacobian_evaluations;
      const size_t number_hessian_evaluations;
      const size_t number_subproblems_solved;

      void print(bool print_primal_dual_solution) const;
   };
} // namespace

#endif // UNO_RESULT_H