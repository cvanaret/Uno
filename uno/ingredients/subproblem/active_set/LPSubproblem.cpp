// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPSubproblem.hpp"
#include "solvers/LP/LPSolverFactory.hpp"

LPSubproblem::LPSubproblem(size_t max_number_variables, size_t max_number_constraints, const Options& options) :
      ActiveSetSubproblem(max_number_variables, max_number_constraints),
      solver(LPSolverFactory::create(max_number_variables, max_number_constraints, options.get_string("LP_solver"), options)) {
}

void LPSubproblem::evaluate_functions(const NonlinearProblem& problem, Iterate& current_iterate) {
   // objective gradient
   problem.evaluate_objective_gradient(current_iterate, this->objective_gradient);

   // constraints
   problem.evaluate_constraints(current_iterate, this->constraints);

   // constraint Jacobian
   problem.evaluate_constraint_jacobian(current_iterate, this->constraint_jacobian);
}

Direction LPSubproblem::solve(Statistics& /*statistics*/, const NonlinearProblem& problem, Iterate& current_iterate) {
   // evaluate the functions at the current iterate
   this->evaluate_functions(problem, current_iterate);

   // bounds of the variable displacements
   this->set_variable_displacement_bounds(problem, current_iterate);

   // bounds of the linearized constraints
   this->set_linearized_constraint_bounds(problem, current_iterate.original_evaluations.constraints);

   return this->solve_LP(problem, current_iterate);
}

Direction LPSubproblem::compute_second_order_correction(const NonlinearProblem& problem, Iterate& trial_iterate) {
   // TODO warm start
   DEBUG << "\nEntered SOC computation\n";
   // shift the RHS with the values of the constraints at the trial iterate
   ActiveSetSubproblem::shift_linearized_constraint_bounds(problem, trial_iterate.original_evaluations.constraints);
   return this->solve_LP(problem, trial_iterate);
}

Direction LPSubproblem::solve_LP(const NonlinearProblem& problem, Iterate& iterate) {
   Direction direction = this->solver->solve_LP(problem.number_variables, problem.number_constraints, this->variable_displacement_bounds,
         this->linearized_constraint_bounds, this->objective_gradient, this->constraint_jacobian, this->initial_point);
   Subproblem::check_unboundedness(direction);
   ActiveSetSubproblem::compute_dual_displacements(problem, iterate, direction);
   this->number_subproblems_solved++;
   return direction;
}

PredictedOptimalityReductionModel LPSubproblem::generate_predicted_optimality_reduction_model(const NonlinearProblem& /*problem*/,
      const Direction& direction) const {
   return PredictedOptimalityReductionModel(-direction.objective, [&]() { // capture "direction" by reference
      // return a function of the step length that cheaply assembles the predicted reduction
      return [=](double step_length) { // capture the expensive quantities by value
         return -step_length * direction.objective;
      };
   });
}

/*
OptimalityMeasureModel LPSubproblem::generate_optimality_measure_model(const ReformulatedProblem& problem, const Direction& direction) const {
   return OptimalityMeasureModel(-direction.objective, [&]() { // capture "direction" by reference
      // return a function of the step length that cheaply assembles the predicted reduction
      return [=](double step_length) { // capture the expensive quantities by value
         return -step_length * direction.objective;
      };
   });
}
*/

size_t LPSubproblem::get_hessian_evaluation_count() const {
   // no second order evaluation is used
   return 0;
}