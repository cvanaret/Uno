#include "LPSubproblem.hpp"
#include "solvers/QP/LPSolverFactory.hpp"

LPSubproblem::LPSubproblem(const NonlinearProblem& problem, size_t max_number_variables, const Options& options) :
      ActiveSetSubproblem(max_number_variables, problem.number_constraints, NO_SOC, false, norm_from_string(options.at("residual_norm"))),
      solver(LPSolverFactory::create(max_number_variables, problem.number_constraints, options.at("LP_solver"))),
      objective_gradient(max_number_variables),
      constraints(problem.number_constraints),
      constraint_jacobian(problem.number_constraints) {
   for (auto& constraint_gradient: this->constraint_jacobian) {
      constraint_gradient.reserve(max_number_variables);
   }
   // register the variables bounds
   for (size_t i = 0; i < problem.number_variables; i++) {
      this->variable_bounds[i] = {problem.get_variable_lower_bound(i), problem.get_variable_upper_bound(i)};
   }
}

void LPSubproblem::evaluate_problem(const NonlinearProblem& problem, Iterate& current_iterate) {
   // objective gradient
   problem.evaluate_objective_gradient(current_iterate, this->objective_gradient);

   // constraints
   problem.evaluate_constraints(current_iterate, this->constraints);

   // constraint Jacobian
   problem.evaluate_constraint_jacobian(current_iterate, this->constraint_jacobian);
}

Direction LPSubproblem::solve(Statistics& /*statistics*/, const NonlinearProblem& problem, Iterate& current_iterate) {
   // evaluate the functions at the current iterate
   this->evaluate_problem(problem, current_iterate);

   // bounds of the variable displacements
   this->set_variable_displacement_bounds(problem, current_iterate);

   // bounds of the linearized constraints
   this->set_linearized_constraint_bounds(problem, current_iterate.original_evaluations.constraints);

   // solve the LP
   Direction direction = this->solver->solve_LP(problem.number_variables, problem.number_constraints, this->variable_displacement_bounds,
         this->linearized_constraint_bounds, this->objective_gradient, this->constraint_jacobian, this->initial_point);
   ActiveSetSubproblem::compute_dual_displacements(problem, current_iterate, direction);
   this->number_subproblems_solved++;
   return direction;
}

PredictedReductionModel LPSubproblem::generate_predicted_reduction_model(const NonlinearProblem& /*problem*/, const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() { // capture direction by reference
      // return a function of the step length that cheaply assembles the predicted reduction
      return [=](double step_length) { // capture the expensive quantities by value
         return -step_length * direction.objective;
      };
   });
}

size_t LPSubproblem::get_hessian_evaluation_count() const {
   // no second order evaluation is used
   return 0;
}

double LPSubproblem::get_proximal_coefficient() const {
   return 0.;
}