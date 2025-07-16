// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "UnoSolverWrapper.hpp"

#include "PythonModel.hpp"
#include "linear_algebra/COOFormat.hpp"
#include "linear_algebra/SparseSymmetricMatrix.hpp"

namespace uno {
   UnoSolverWrapper::UnoSolverWrapper(bool constrained_model, const Options& options):
         uno_solver(constrained_model, options) {
   }

   void UnoSolverWrapper::solve(size_t number_variables, size_t number_constraints, const objective_function_type& evaluate_objective,
         const constraint_functions_type& evaluate_constraints, const objective_gradient_type& evaluate_objective_gradient,
         const jacobian_type& evaluate_jacobian, const lagrangian_hessian_type& evaluate_lagrangian_hessian,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds,
         const Options& options) {
      std::cout << "Congrats, you just called Uno with (n, m) = (" << number_variables << ", " << number_constraints << ")\n";

      /*
      // hs015 primal-dual iterate
      Vector<double> x{-2., 1.};
      Vector<double> y{0., 0.};
      const double objective_multiplier = 1.;

      // data structures
      Vector<double> constraints(number_constraints);
      SparseVector<double> objective_gradient(number_variables);
      SparseSymmetricMatrix<COOFormat<size_t, double>> jacobian(number_variables, 4, 0);
      SparseSymmetricMatrix<COOFormat<size_t, double>> hessian(number_variables, 3, 0);

      // evaluations
      const double objective = evaluate_objective(wrap_pointer(&x));
      evaluate_constraints(wrap_pointer(&x), wrap_pointer(&constraints));
      evaluate_objective_gradient(wrap_pointer(&x), wrap_pointer(&objective_gradient));
      evaluate_jacobian(wrap_pointer(&x), wrap_pointer<SymmetricMatrix<size_t, double>>(&jacobian));
      evaluate_lagrangian_hessian(wrap_pointer(&x), objective_multiplier, wrap_pointer(&y),
         wrap_pointer<SymmetricMatrix<size_t, double>>(&hessian));

      std::cout << "Evaluations at x = " << x << "\n";
      std::cout << "Objective value: " << objective << '\n';
      std::cout << "Constraint values: " << constraints << '\n';
      std::cout << "Objective gradient: " << objective_gradient;
      std::cout << "Jacobian: " << jacobian;
      std::cout << "Hessian: " << hessian;
      */

      const PythonModel model("python_model", number_variables, number_constraints,
         1., evaluate_objective, evaluate_constraints, evaluate_objective_gradient, evaluate_jacobian,
         evaluate_lagrangian_hessian, variables_lower_bounds, variables_upper_bounds, constraints_lower_bounds,
         constraints_upper_bounds);

      // initialize initial primal and dual points
      Iterate initial_iterate(number_variables, number_constraints);
      model.initial_primal_point(initial_iterate.primals);
      model.project_onto_variable_bounds(initial_iterate.primals);
      model.initial_dual_point(initial_iterate.multipliers.constraints);
      initial_iterate.feasibility_multipliers.reset();

      // solve the instance
      const Result result = this->uno_solver.solve(model, initial_iterate, options);
      if (result.optimization_status == OptimizationStatus::SUCCESS) {
         // check result.solution.status
      }
      else {
         // ...
      }
   }
} // namespace