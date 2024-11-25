// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RESULT_H
#define UNO_RESULT_H

#include "Iterate.hpp"
#include "OptimizationStatus.hpp"

namespace uno {
   struct Result {
      Result() = delete;

      OptimizationStatus optimization_status;
      Iterate solution;
      size_t number_variables;
      size_t number_constraints;
      size_t iteration;
      double cpu_time;
      size_t objective_evaluations;
      size_t constraint_evaluations;
      size_t objective_gradient_evaluations;
      size_t jacobian_evaluations;
      size_t hessian_evaluations;
      size_t number_subproblems_solved;

      void print(bool print_primal_dual_solution) const;
   };
} // namespace

#endif // UNO_RESULT_H
