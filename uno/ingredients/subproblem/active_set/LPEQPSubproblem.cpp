// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPEQPSubproblem.hpp"
#include "solvers/QP/QPSolverFactory.hpp"
#include "tools/Infinity.hpp"

// Question: Can we switch Hessian off easily in QP solve?
// Easy for BQPD: simply return 0 in gdotx (not super-efficient)
// ===========================================================
LPEQPSubproblem::LPEQPSubproblem(size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros, const Options& options) :
      ActiveSetSubproblem(max_number_variables, max_number_constraints),
      use_regularization(true),
      // if no trust region is used, the problem should be convexified to guarantee boundedness + descent direction
      hessian_model(HessianModelFactory::create(options.get_string("hessian_model"), max_number_variables,
            max_number_hessian_nonzeros + max_number_variables, options.get_string("mechanism") != "TR", options)),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory::create(options.get_string("QP_solver"), max_number_variables, max_number_constraints,
            hessian_model->hessian->capacity, true, options)),
      statistics_regularization_column_order(options.get_int("statistics_regularization_column_order")) {
}

void LPEQPSubproblem::initialize(Statistics& statistics, const NonlinearProblem& /*problem*/, Iterate& /*first_iterate*/) {
   if (this->use_regularization) {
      statistics.add_column("regularization", Statistics::double_width, this->statistics_regularization_column_order);
   }
}

void LPEQPSubproblem::set_variable_EQP_bounds(const NonlinearProblem& problem, const Iterate& current_iterate, Direction& LP_direction) {
   // initialize all bounds to {-INF<double>,+INF<double>}
   for (size_t i: Range(problem.number_variables)) {
     this->variable_displacement_bounds[i] = {-INF<double>,+INF<double>};
   }
   // set active lower bounds as equality constraints 
   for (size_t i: LP_direction.active_set.bounds.at_lower_bound) {
      const double lb = this->variable_bounds[i].lb - current_iterate.primals[i];
      if (lb >= -(this->trust_region_radius)) {
	this->variable_displacement_bounds[i] = {lb, lb};
      }
   }
   // set active lower bounds as equality constraints 
   for (size_t i: LP_direction.active_set.bounds.at_upper_bound) {
      const double ub = this->variable_bounds[i].ub - current_iterate.primals[i];
      if (ub <= this->trust_region_radius) {
	this->variable_displacement_bounds[i] = {ub, ub};
      }
   }
}

void LPEQPSubproblem::set_linearized_EQP_bounds(const NonlinearProblem& problem, const std::vector<double>& current_constraints, Direction& LP_direction) {
   // initialize all constraint bounds to {-INF<double>,+INF<double>}
   for (size_t j: Range(problem.number_constraints)) {
     this->linearized_constraint_bounds[j] = {-INF<double>,+INF<double>};
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

void LPEQPSubproblem::evaluate_functions(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate) {
   // Hessian
   this->hessian_model->evaluate(statistics, problem, current_iterate.primals, current_iterate.multipliers.constraints);
   // objective gradient, constraints and constraint Jacobian
   problem.evaluate_objective_gradient(current_iterate, this->evaluations.objective_gradient);
   problem.evaluate_constraints(current_iterate, this->evaluations.constraints);
   problem.evaluate_constraint_jacobian(current_iterate, this->evaluations.constraint_jacobian);
}

Direction LPEQPSubproblem::solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate) {
    
   // evaluate the functions at the current iterate
   this->evaluate_functions(statistics, problem, current_iterate);

   // bounds of the variable displacements
   this->set_variable_bounds(problem, current_iterate);
   this->set_variable_displacement_bounds(problem, current_iterate);

   // bounds of the linearized constraints
   this->set_linearized_constraint_bounds(problem, this->evaluations.constraints);

   // solve LP subproblem
   Direction LP_direction = this->solve_LP(problem, current_iterate);

   DEBUG << "d^*(LP) = ";
   print_vector(DEBUG, LP_direction.primals, 0, LP_direction.number_variables);
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

   // set-up EQP subproblem: set inactive bounds +/- INF and others as equations
   this->evaluate_functions(statistics, problem, current_iterate);
   this->set_variable_EQP_bounds(problem, current_iterate, LP_direction);
   this->set_linearized_EQP_bounds(problem, this->evaluations.constraints, LP_direction);

   // return EQP solution
   return this->solve_QP(problem, current_iterate);
}

Direction LPEQPSubproblem::compute_second_order_correction(const NonlinearProblem& /*problem*/, Iterate& /*trial_iterate*/) {
   // TODO warm start
   DEBUG << "\nEntered SOC computation\n";
   assert(false && "Not implemented yet");
}

Direction LPEQPSubproblem::solve_QP(const NonlinearProblem& problem, Iterate& iterate) {
   Direction direction = this->solver->solve_QP(problem.number_variables, problem.number_constraints, this->variable_displacement_bounds,
         this->linearized_constraint_bounds, this->evaluations.objective_gradient, this->evaluations.constraint_jacobian,
         *this->hessian_model->hessian, this->initial_point);
   Subproblem::check_unboundedness(direction);
   ActiveSetSubproblem::compute_dual_displacements(problem, iterate, direction);
   this->number_subproblems_solved++;
   return direction;
}

Direction LPEQPSubproblem::solve_LP(const NonlinearProblem& problem, Iterate& iterate) {
   Direction direction = this->solver->solve_LP(problem.number_variables, problem.number_constraints, this->variable_displacement_bounds,
         this->linearized_constraint_bounds, this->evaluations.objective_gradient, this->evaluations.constraint_jacobian,
         this->initial_point);
   Subproblem::check_unboundedness(direction);
   ActiveSetSubproblem::compute_dual_displacements(problem, iterate, direction);
   this->number_subproblems_solved++;
   return direction;
}

size_t LPEQPSubproblem::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}
