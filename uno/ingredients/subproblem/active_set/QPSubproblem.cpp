// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "QPSubproblem.hpp"
#include "solvers/QP/QPSolverFactory.hpp"

QPSubproblem::QPSubproblem(Statistics& statistics, size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros,
         const Options& options) :
      ActiveSetSubproblem(max_number_variables, max_number_constraints),
      use_regularization(options.get_string("globalization_mechanism") != "TR"),
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

void QPSubproblem::generate_initial_iterate(const NonlinearProblem& /*problem*/, Iterate& /*initial_iterate*/) {
}

void QPSubproblem::evaluate_functions(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate) {
   // Hessian
   this->hessian_model->evaluate(statistics, problem, current_iterate.primals, current_iterate.multipliers.constraints);
   // objective gradient, constraints and constraint Jacobian
   problem.evaluate_objective_gradient(current_iterate, this->evaluations.objective_gradient);
   problem.evaluate_constraints(current_iterate, this->evaluations.constraints);
   problem.evaluate_constraint_jacobian(current_iterate, this->evaluations.constraint_jacobian);
}

Direction QPSubproblem::solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate, bool evaluate_functions) {
   if (evaluate_functions) {
      // evaluate the functions at the current iterate
      this->evaluate_functions(statistics, problem, current_iterate);
   }

   // bounds of the variable displacements
   this->set_variable_bounds(problem, current_iterate);
   this->set_variable_displacement_bounds(problem, current_iterate);

   // bounds of the linearized constraints
   this->set_linearized_constraint_bounds(problem, this->evaluations.constraints);

   return this->solve_QP(problem, current_iterate);
}

Direction QPSubproblem::solve_QP(const NonlinearProblem& problem, Iterate& iterate) {
   Direction direction = this->solver->solve_QP(problem.number_variables, problem.number_constraints, this->variable_displacement_bounds,
         this->linearized_constraint_bounds, this->evaluations.objective_gradient, this->evaluations.constraint_jacobian,
         *this->hessian_model->hessian, this->initial_point);
   Subproblem::check_unboundedness(direction);
   ActiveSetSubproblem::compute_dual_displacements(problem, iterate, direction);
   this->number_subproblems_solved++;
   return direction;
}

size_t QPSubproblem::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}