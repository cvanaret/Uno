// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "QPSubproblem.hpp"
#include "solvers/QP/QPSolverFactory.hpp"

QPSubproblem::QPSubproblem(Statistics& statistics, size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros,
         const Options& options) :
      ActiveSetSubproblem(max_number_variables, max_number_constraints),
      use_regularization(options.get_string("globalization_mechanism") != "TR" || options.get_bool("convexify_QP")),
      // if no trust region is used, the problem should be convexified to guarantee boundedness
      hessian_model(HessianModelFactory::create(options.get_string("hessian_model"), max_number_variables,
            max_number_hessian_nonzeros + max_number_variables, this->use_regularization, options)),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory::create(options.get_string("QP_solver"), max_number_variables, max_number_constraints,
            hessian_model->hessian->capacity, true, options)) {
   if (this->use_regularization) {
      statistics.add_column("regularization", Statistics::double_width, options.get_int("statistics_regularization_column_order"));
   }
}

void QPSubproblem::evaluate_functions(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
      const WarmstartInformation& warmstart_information) {
   // Lagrangian Hessian
   if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
      this->hessian_model->evaluate(statistics, problem, current_iterate.primals, current_iterate.multipliers.constraints);
   }
   // objective gradient, constraints and constraint Jacobian
   if (warmstart_information.objective_changed) {
      problem.evaluate_objective_gradient(current_iterate, this->evaluations.objective_gradient);
   }
   if (warmstart_information.constraints_changed) {
      problem.evaluate_constraints(current_iterate, this->evaluations.constraints);
      problem.evaluate_constraint_jacobian(current_iterate, this->evaluations.constraint_jacobian);
   }
}

Direction QPSubproblem::solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
      const WarmstartInformation& warmstart_information) {
   // evaluate the functions at the current iterate
   this->evaluate_functions(statistics, problem, current_iterate, warmstart_information);

   // set bounds of the variable displacements
   if (warmstart_information.variable_bounds_changed) {
      this->set_direction_bounds(problem, current_iterate);
   }

   // set bounds of the linearized constraints
   if (warmstart_information.constraint_bounds_changed) {
      this->set_linearized_constraint_bounds(problem, this->evaluations.constraints);
   }

   // solve the QP
   Direction direction = this->solver->solve_QP(problem.number_variables, problem.number_constraints, this->direction_bounds,
         this->linearized_constraint_bounds, this->evaluations.objective_gradient, this->evaluations.constraint_jacobian,
         *this->hessian_model->hessian, this->initial_point, warmstart_information);
   ActiveSetSubproblem::compute_dual_displacements(problem, current_iterate, direction);
   this->number_subproblems_solved++;
   // reset the initial point
   initialize_vector(this->initial_point, 0.);
   return direction;
}

std::function<double(double)> QPSubproblem::compute_predicted_optimality_reduction_model(const NonlinearProblem& problem,
      const Iterate& current_iterate, const Direction& direction, double step_length) const {
   return problem.compute_predicted_optimality_reduction_model(current_iterate, direction, step_length, *this->hessian_model->hessian);
}

size_t QPSubproblem::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}