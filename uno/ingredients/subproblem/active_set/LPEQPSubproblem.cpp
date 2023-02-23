// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Question: Do we need to declare an LP solver too? Possible?
// ===========================================================
#include "LPEQPSubproblem.hpp"
#include "solvers/QP/QPSolverFactory.hpp"
#include "solvers/LP/LPSolverFactory.hpp"

// Question: Can we switch Hessian off easily in QP solve?
// Easy for BQPD: simply return 0 in gdotx (not super-efficient)
// ===========================================================
LPEQPSubproblem::LPEQPSubproblem(size_t max_number_variables, size_t max_number_constraints, size_t max_number_hessian_nonzeros, const Options& options) :
      ActiveSetSubproblem(max_number_variables, max_number_constraints),
      // if no trust region is used, the problem should be convexified to guarantee boundedness + descent direction
      hessian_model(HessianModelFactory::create(options.get_string("hessian_model"), max_number_variables,
            max_number_hessian_nonzeros + max_number_variables, options.get_string("mechanism") != "TR", options)),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory::create(options.get_string("QP_solver"), max_number_variables, max_number_constraints,
            hessian_model->hessian->capacity, true, options)) {
}

void LPEQPSubproblem::evaluate_functionsEQP(const NonlinearProblem& problem, Iterate& current_iterate) {
   // Hessian
   this->hessian_model->evaluate(problem, current_iterate.primals, current_iterate.multipliers.constraints);
   // objective gradient, constraints and constraint Jacobian
   problem.evaluate_objective_gradient(current_iterate, this->evaluations.objective_gradient);
   problem.evaluate_constraints(current_iterate, this->evaluations.constraints);
   problem.evaluate_constraint_jacobian(current_iterate, this->evaluations.constraint_jacobian);
}

void LPEQPSubproblem::evaluate_functionsLP(const NonlinearProblem& problem, Iterate& current_iterate) {
   // objective gradient, constraints and constraint Jacobian
   // objective gradient, constraints and constraint Jacobian
   problem.evaluate_objective_gradient(current_iterate, this->evaluations.objective_gradient);
   problem.evaluate_constraints(current_iterate, this->evaluations.constraints);
   problem.evaluate_constraint_jacobian(current_iterate, this->evaluations.constraint_jacobian);
}

Direction LPEQPSubproblem::solve(Statistics& /*statistics*/, const NonlinearProblem& problem, Iterate& current_iterate) {
    
  // evaluate the functions at the current iterate
   this->evaluate_functionsLP(problem, current_iterate);

   // bounds of the variable displacements
   this->set_variable_bounds(problem, current_iterate);
   this->set_variable_displacement_bounds(problem, current_iterate);

   // bounds of the linearized constraints
   this->set_linearized_constraint_bounds(problem, this->evaluations.constraints);

   // solve LP subproblem
   Direction LP_direction = this->solve_LP(problem, current_iterate);
   
   // set-up EQP subproblem:
   // (1) set all bnds=+/-\infty;
   // (2) set active constraints as equations 

   

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

Direction LPSubproblem::solve_LP(const NonlinearProblem& problem, Iterate& iterate) {
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
