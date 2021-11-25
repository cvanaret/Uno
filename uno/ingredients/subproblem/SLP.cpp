#include "SLP.hpp"
#include "solvers/QP/LPSolverFactory.hpp"

SLP::SLP(const Problem& problem, size_t max_number_variables, size_t number_constraints, const std::string& LP_solver_name) :
      Subproblem(problem.number_variables, max_number_variables, number_constraints, NO_SOC, false),
      solver(LPSolverFactory::create(max_number_variables, number_constraints, LP_solver_name)),
      initial_point(max_number_variables) {
   // register the original constraints bounds
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->constraints_bounds[j] = problem.constraint_bounds[j];
   }
}

void SLP::create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);
   /* compute first- and second-order information */
   problem.evaluate_constraints(current_iterate.x, current_iterate.constraints);
   for (auto& row: this->constraints_jacobian) {
      row.clear();
   }
   problem.evaluate_constraint_jacobian(current_iterate.x, this->constraints_jacobian);

   this->build_objective_model(problem, current_iterate, objective_multiplier);

   /* bounds of the variables */
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->set_constraints_bounds(problem, current_iterate.constraints);

   /* set the initial point */
   clear(this->initial_point);
}

void SLP::build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // objective gradient
   this->set_scaled_objective_gradient(problem, current_iterate, objective_multiplier);

   // initial point
   clear(this->initial_point);
}

void SLP::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

Direction SLP::solve(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) {
   /* solve the LP */
   Direction direction = this->solver->solve_LP(this->variables_bounds, this->constraints_bounds, this->objective_gradient,
         this->constraints_jacobian, this->initial_point);
   this->number_subproblems_solved++;

   // compute dual displacements (SLP methods usually compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
   return direction;
}

PredictedReductionModel SLP::generate_predicted_reduction_model(const Problem& /*problem*/, const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() { // capture direction by reference
      // return a function of the step length that cheaply assembles the predicted reduction
      return [=](double step_length) { // capture the expensive quantities by value
         return -step_length * direction.objective;
      };
   });
}

size_t SLP::get_hessian_evaluation_count() const {
   // no second order evaluation is used
   return 0;
}