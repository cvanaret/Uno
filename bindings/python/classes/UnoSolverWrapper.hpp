// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOSOLVERWRAPPER_H
#define UNO_UNOSOLVERWRAPPER_H


#include "Uno.hpp"
#include "PythonTypes.hpp"

namespace uno {
   // forward declarations
   class Options;

   class UnoSolverWrapper {
   public:
      UnoSolverWrapper(bool constrained_model, const Options& options);

      void solve(size_t number_variables, size_t number_constraints, const objective_function_type& evaluate_objective,
         const constraint_functions_type& evaluate_constraints, const objective_gradient_type& evaluate_objective_gradient,
         const jacobian_type& evaluate_jacobian, const lagrangian_hessian_type& evaluate_lagrangian_hessian,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds,
         const Options& options);

   protected:
      Uno uno_solver;
   };
} // namespace

#endif // UNO_UNOSOLVERWRAPPER_H