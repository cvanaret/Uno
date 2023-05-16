// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPEQPSubproblem.hpp"
#include "solvers/QP/QPSolverFactory.hpp"
#include "tools/Infinity.hpp"

// Question: Can we switch Hessian off easily in QP solve?
// Easy for BQPD: simply return 0 in gdotx (not super-efficient)
// ===========================================================
LPEQPSubproblem::LPEQPSubproblem(Statistics& statistics, size_t max_number_variables, size_t max_number_constraints,
         size_t max_number_hessian_nonzeros, const Options& options) :
      ActiveSetSubproblem(max_number_variables, max_number_constraints),
      use_regularization(true),
      // if no trust region is used, the problem should be convexified to guarantee boundedness + descent direction
      hessian_model(HessianModelFactory::create(options.get_string("hessian_model"), max_number_variables,
            max_number_hessian_nonzeros + max_number_variables, this->use_regularization, options)),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory::create(options.get_string("QP_solver"), max_number_variables, max_number_constraints,
            hessian_model->hessian->capacity, true, options)) {
   if (this->use_regularization) {
      statistics.add_column("regularization", Statistics::double_width, options.get_int("statistics_regularization_column_order"));
   }
}

void LPEQPSubproblem::generate_initial_iterate(const NonlinearProblem& /*problem*/, Iterate& /*first_iterate*/) {
}

void LPEQPSubproblem::evaluate_functions(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
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

Direction LPEQPSubproblem::solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
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

   // solve LP subproblem
   Direction LP_direction = this->solve_LP(problem, current_iterate, warmstart_information);

   DEBUG << "d^*(LP) = "; print_vector(DEBUG, LP_direction.primals, 0, LP_direction.number_variables);
   DEBUG << "bound constraints active at lower bound =";
   for (size_t i: LP_direction.active_set.bounds.at_lower_bound) {
      DEBUG << " x" << i;
   }
   DEBUG << '\n';
   DEBUG << "bound constraints active at upper bound =";
   for (size_t i: LP_direction.active_set.bounds.at_upper_bound) {
      DEBUG << " x" << i;
   }
   DEBUG << '\n';
   DEBUG << "constraints at lower bound =";
   for (size_t j: LP_direction.active_set.constraints.at_lower_bound) {
      DEBUG << " c" << j;
   }
   DEBUG << '\n';
   DEBUG << "constraints at upper bound =";
   for (size_t j: LP_direction.active_set.constraints.at_upper_bound) {
      DEBUG << " c" << j;
   }
   DEBUG << '\n';

   // set up EQP subproblem: set inactive bounds +/- INF and others as equations
   this->set_variable_EQP_bounds(problem, current_iterate, LP_direction);
   this->set_linearized_EQP_bounds(problem, this->evaluations.constraints, LP_direction);

   // return EQP solution
   return this->solve_QP(problem, current_iterate, warmstart_information);
}

Direction LPEQPSubproblem::solve_LP(const NonlinearProblem& problem, Iterate& iterate, const WarmstartInformation& warmstart_information) {
   Direction direction = this->solver->solve_LP(problem.number_variables, problem.number_constraints, this->direction_bounds,
         this->linearized_constraint_bounds, this->evaluations.objective_gradient, this->evaluations.constraint_jacobian,
         this->initial_point, warmstart_information);
   ActiveSetSubproblem::compute_dual_displacements(problem, iterate, direction);
   this->number_subproblems_solved++;
   return direction;
}

Direction LPEQPSubproblem::solve_QP(const NonlinearProblem& problem, Iterate& iterate, const WarmstartInformation& warmstart_information) {
   Direction direction = this->solver->solve_QP(problem.number_variables, problem.number_constraints, this->direction_bounds,
         this->linearized_constraint_bounds, this->evaluations.objective_gradient, this->evaluations.constraint_jacobian,
         *this->hessian_model->hessian, this->initial_point, warmstart_information);
   ActiveSetSubproblem::compute_dual_displacements(problem, iterate, direction);
   this->number_subproblems_solved++;
   return direction;
}

void LPEQPSubproblem::set_variable_EQP_bounds(const NonlinearProblem& problem, const Iterate& current_iterate, Direction& LP_direction) {
   // initialize all bounds to {-INF<double>,+INF<double>}
   for (size_t i: Range(problem.number_variables)) {
      this->direction_bounds[i] = {-INF<double>, +INF<double>};
   }
   // set active lower bounds as equality constraints
   for (size_t i: LP_direction.active_set.bounds.at_lower_bound) {
      const double lb = problem.get_variable_lower_bound(i) - current_iterate.primals[i];
      if (lb >= -(this->trust_region_radius)) {
         this->direction_bounds[i] = {lb, lb};
      }
   }
   // set active lower bounds as equality constraints
   for (size_t i: LP_direction.active_set.bounds.at_upper_bound) {
      const double ub = problem.get_variable_upper_bound(i) - current_iterate.primals[i];
      if (ub <= this->trust_region_radius) {
         this->direction_bounds[i] = {ub, ub};
      }
   }
}

void LPEQPSubproblem::set_linearized_EQP_bounds(const NonlinearProblem& problem, const std::vector<double>& current_constraints, Direction& LP_direction) {
   // initialize all constraint bounds to {-INF<double>,+INF<double>}
   for (size_t j: Range(problem.number_constraints)) {
      this->linearized_constraint_bounds[j] = {-INF<double>, +INF<double>};
   }
   // set active lower constraint bounds as equality constraints
   for (size_t j: LP_direction.active_set.constraints.at_lower_bound) {
      const double lb = problem.get_constraint_lower_bound(j) - current_constraints[j];
      this->linearized_constraint_bounds[j] = {lb, lb};
   }
   // set active upper constraint bounds as equality constraints
   for (size_t j: LP_direction.active_set.constraints.at_upper_bound) {
      const double ub = problem.get_constraint_upper_bound(j) - current_constraints[j];
      this->linearized_constraint_bounds[j] = {ub, ub};
   }
}

size_t LPEQPSubproblem::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}
