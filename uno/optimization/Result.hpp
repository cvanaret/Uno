// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RESULT_H
#define UNO_RESULT_H

#include "Iterate.hpp"
#include "TerminationStatus.hpp"

struct Result {
   Result() = delete;

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

#endif // UNO_RESULT_H