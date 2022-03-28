#include "LPSubproblem.hpp"
#include "solvers/QP/LPSolverFactory.hpp"

LPSubproblem::LPSubproblem(const Problem& problem, size_t max_number_variables, const Options& options) :
      ActiveSetSubproblem(max_number_variables, problem.number_constraints, NO_SOC, false, norm_from_string(options.at("residual_norm"))),
      solver(LPSolverFactory::create(max_number_variables, problem.number_constraints, options.at("LP_solver"))) {
}

void LPSubproblem::build_objective_model(const Problem& problem, Iterate& current_iterate, double /*objective_multiplier*/) {
   // objective gradient
   current_iterate.evaluate_objective_gradient(problem);
}

void LPSubproblem::build_constraint_model(const Problem& problem, Iterate& current_iterate) {
   // compute first- and second-order information
   current_iterate.evaluate_constraints(problem);
   current_iterate.evaluate_constraint_jacobian(problem);
}

Direction LPSubproblem::solve(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) {
   // bounds of the variable displacements
   this->set_variable_displacement_bounds(problem, current_iterate);

   // bounds of the linearized constraints
   this->set_linearized_constraint_bounds(problem, current_iterate.problem_evaluations.constraints);

   // solve the LP
   Direction direction = this->solver->solve_LP(problem.number_variables, problem.number_constraints, this->variable_displacement_bounds,
         this->linearized_constraint_bounds, current_iterate.subproblem_evaluations.objective_gradient,
         current_iterate.problem_evaluations.constraint_jacobian, this->initial_point);
   ActiveSetSubproblem::compute_dual_displacements(problem, current_iterate, direction);
   this->number_subproblems_solved++;
   return direction;
}

PredictedReductionModel LPSubproblem::generate_predicted_reduction_model(const Problem& /*problem*/, const Iterate& /*current_iterate*/,
      const Direction& direction) const {
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