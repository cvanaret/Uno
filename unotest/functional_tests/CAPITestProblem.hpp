// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CAPITESTPROBLEM_H
#define UNO_CAPITESTPROBLEM_H

// problem hs015 built through the C API, for the tests of the C API

#include <cmath>
#include "../../interfaces/C/Uno_C_API.h"
#include "tools/Infinity.hpp"

namespace uno::capi_test {
   constexpr uno_int number_variables = 2;
   constexpr uno_int number_constraints = 2;

   inline uno_int objective_function(uno_int /*number_variables*/, const double* x, double* objective_value, void* /*user_data*/) {
      *objective_value = 100. * std::pow(x[1] - x[0]*x[0], 2.) + std::pow(1. - x[0], 2.);
      return 0;
   }

   inline uno_int objective_gradient(uno_int /*number_variables*/, const double* x, double* gradient, void* /*user_data*/) {
      gradient[0] = 400.*std::pow(x[0], 3.) - 400.*x[0]*x[1] + 2.*x[0] - 2.;
      gradient[1] = 200.*(x[1] - x[0]*x[0]);
      return 0;
   }

   inline uno_int constraint_functions(uno_int /*number_variables*/, uno_int /*number_constraints*/, const double* x,
         double* constraint_values, void* /*user_data*/) {
      constraint_values[0] = x[0]*x[1];
      constraint_values[1] = x[0] + x[1]*x[1];
      return 0;
   }

   inline uno_int jacobian(uno_int /*number_variables*/, uno_int /*number_jacobian_nonzeros*/, const double* x,
         double* jacobian_values, void* /*user_data*/) {
      jacobian_values[0] = x[1];
      jacobian_values[1] = 1.;
      jacobian_values[2] = x[0];
      jacobian_values[3] = 2.*x[1];
      return 0;
   }

   inline uno_int lagrangian_hessian(uno_int /*number_variables*/, uno_int /*number_constraints*/, uno_int /*number_hessian_nonzeros*/,
         const double* x, double objective_multiplier, const double* multipliers, double* hessian_values, void* /*user_data*/) {
      hessian_values[0] = objective_multiplier*(1200.*x[0]*x[0] - 400.*x[1] + 2.);
      hessian_values[1] = -400.*objective_multiplier*x[0] - multipliers[0];
      hessian_values[2] = 200.*objective_multiplier - 2.*multipliers[1];
      return 0;
   }

   // a preset whose subproblem solver was compiled in; nullptr if none
   inline const char* available_preset() {
#if defined(HAS_BQPD)
      return "filtersqp";
#elif defined(HAS_HSL) || defined(HAS_MA57) || defined(HAS_MA27) || defined(HAS_MUMPS) || defined(HAS_SPRAL)
      return "ipopt";
#elif defined(HAS_HIGHS)
      return "filterslp";
#else
      return nullptr;
#endif
   }

   // creates the C API model of hs015 (to be destroyed with uno_destroy_model)
   inline uno_model* create_hs015_model() {
      double variables_lower_bounds[] = {-Inf, -Inf};
      double variables_upper_bounds[] = {0.5, Inf};
      double constraints_lower_bounds[] = {1., 0.};
      double constraints_upper_bounds[] = {Inf, Inf};
      uno_int jacobian_row_indices[] = {0, 1, 0, 1};
      uno_int jacobian_column_indices[] = {0, 0, 1, 1};
      uno_int hessian_row_indices[] = {0, 1, 1};
      uno_int hessian_column_indices[] = {0, 0, 1};
      double x0[] = {-2., 1.};

      uno_model* model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds,
         variables_upper_bounds, UNO_ZERO_BASED_INDEXING);
      uno_set_objective(model, UNO_MINIMIZE, objective_function, objective_gradient);
      uno_set_constraints(model, number_constraints, constraint_functions, constraints_lower_bounds,
         constraints_upper_bounds, 4, jacobian_row_indices, jacobian_column_indices, jacobian);
      uno_set_lagrangian_hessian(model, 3, UNO_LOWER_TRIANGLE, hessian_row_indices, hessian_column_indices,
         lagrangian_hessian);
      uno_set_lagrangian_sign_convention(model, UNO_MULTIPLIER_NEGATIVE);
      uno_set_initial_primal_iterate(model, x0);
      return model;
   }
} // namespace

#endif // UNO_CAPITESTPROBLEM_H
