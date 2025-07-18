// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "UnoSolverWrapper.hpp"
#include "PythonModel.hpp"
#include "model/ModelFactory.hpp"
#include "options/Options.hpp"

namespace uno {
   UnoSolverWrapper::UnoSolverWrapper(bool constrained_model, const Options& options):
         uno_solver(constrained_model, options) {
   }

   void UnoSolverWrapper::solve(size_t number_variables, size_t number_constraints, const objective_function_type& evaluate_objective,
         const constraint_functions_type& evaluate_constraints, const objective_gradient_type& evaluate_objective_gradient,
         const jacobian_type& evaluate_jacobian, const lagrangian_hessian_type& evaluate_lagrangian_hessian,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds,
         const std::vector<double>& primal_initial_point, const std::vector<double>& dual_initial_point,
         const Options& options) {
      // create the model
      std::unique_ptr<Model> python_model = std::make_unique<PythonModel>("python_model", number_variables, number_constraints,
         1., evaluate_objective, evaluate_constraints, evaluate_objective_gradient, evaluate_jacobian,
         evaluate_lagrangian_hessian, variables_lower_bounds, variables_upper_bounds, constraints_lower_bounds,
         constraints_upper_bounds, primal_initial_point, dual_initial_point);
      // reformulate (scale, add slacks, relax the bounds, ...) if necessary
      std::unique_ptr<Model> model = ModelFactory::reformulate(std::move(python_model), options);
      DISCRETE << "Reformulated model " << model->name << '\n' << model->number_variables << " variables, " <<
         model->number_constraints << " constraints (" << model->get_equality_constraints().size() <<
         " equality, " << model->get_inequality_constraints().size() << " inequality)\n";

      Logger::set_logger(options.get_string("logger"));
      DISCRETE << "Original model " << model->name << '\n' << model->number_variables << " variables, " <<
         model->number_constraints << " constraints (" << model->get_equality_constraints().size() <<
         " equality, " << model->get_inequality_constraints().size() << " inequality)\n";

      // initialize initial primal and dual points
      Iterate initial_iterate(number_variables, number_constraints);
      model->initial_primal_point(initial_iterate.primals);
      model->project_onto_variable_bounds(initial_iterate.primals);
      model->initial_dual_point(initial_iterate.multipliers.constraints);
      initial_iterate.feasibility_multipliers.reset();

      // solve the instance
      const Result result = this->uno_solver.solve(*model, initial_iterate, options);
      if (result.optimization_status == OptimizationStatus::SUCCESS) {
         // check result.solution.status
      }
      else {
         // ...
      }
   }
} // namespace