#include "LPSubproblem.hpp"
#include "solvers/QP/LPSolverFactory.hpp"

LPSubproblem::LPSubproblem(const Problem& problem, size_t max_number_variables, const Options& options) :
      Subproblem(problem.number_variables, max_number_variables, problem.number_constraints, false, NO_SOC, false,
            norm_from_string(options.at("residual_norm"))),
      solver(LPSolverFactory::create(max_number_variables, this->number_constraints, options.at("LP_solver"))),
      initial_point(max_number_variables) {
}

void LPSubproblem::build_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier,
      double trust_region_radius) {
   copy_from(this->constraint_multipliers, current_iterate.multipliers.constraints);
   // compute first- and second-order information
   problem.evaluate_constraints(current_iterate.x, current_iterate.constraints);
   for (auto& row: this->constraint_jacobian) {
      row.clear();
   }
   current_iterate.evaluate_constraint_jacobian(problem);
   this->constraint_jacobian = current_iterate.constraint_jacobian;

   this->build_objective_model(problem, current_iterate, objective_multiplier);

   // bounds of the variables
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   // bounds of the linearized constraints
   this->set_constraint_bounds(problem, current_iterate.constraints);

   // set the initial point
   initialize_vector(this->initial_point, 0.);
}

void LPSubproblem::build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // objective gradient
   this->set_scaled_objective_gradient(problem, current_iterate, objective_multiplier);

   // initial point
   initialize_vector(this->initial_point, 0.);
}

void LPSubproblem::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

Direction LPSubproblem::solve(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) {
   // solve the LP
   Direction direction = this->solver->solve_LP(this->variables_bounds, this->constraint_bounds, this->objective_gradient,
         this->constraint_jacobian, this->initial_point);
   this->number_subproblems_solved++;

   // compute dual displacements (SLP methods usually compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
   return direction;
}

PredictedReductionModel LPSubproblem::generate_predicted_reduction_model(const Problem& /*problem*/, const Direction& direction) const {
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