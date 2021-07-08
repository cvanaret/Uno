#include "SLP.hpp"
#include "QPSolverFactory.hpp"

SLP::SLP(size_t number_variables, size_t number_constraints, const std::string& QP_solver_name) :
      Subproblem(number_variables, number_constraints),
      solver(QPSolverFactory::create(QP_solver_name, number_variables, number_constraints, 0, false)),
      initial_point(number_variables) {
}

void SLP::generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);
   /* compute first- and second-order information */
   problem.evaluate_constraints(current_iterate.x, current_iterate.constraints);
   for (auto& row: this->constraints_jacobian) {
      row.clear();
   }
   problem.constraints_jacobian(current_iterate.x, this->constraints_jacobian);

   this->objective_gradient.clear();
   problem.objective_gradient(current_iterate.x, this->objective_gradient);
   this->update_objective_multiplier(problem, current_iterate, objective_multiplier);

   /* bounds of the variables */
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->set_constraints_bounds(problem, current_iterate.constraints);

   /* set the initial point */
   clear(this->initial_point);
}

void SLP::update_objective_multiplier(const Problem& /*problem*/, const Iterate& current_iterate, double objective_multiplier) {
   // scale objective gradient
   if (objective_multiplier == 0.) {
      clear(this->objective_gradient);
   }
   else if (objective_multiplier < 1.) {
      this->objective_gradient = current_iterate.objective_gradient;
      scale(this->objective_gradient, objective_multiplier);
   }
   clear(this->initial_point);
}

void SLP::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

Direction SLP::compute_direction(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) {
   /* solve the LP */
   Direction direction = this->solver->solve_LP(variables_bounds, constraints_bounds, this->objective_gradient,
         this->constraints_jacobian,this->initial_point);
   // compute dual displacements (SQP methods compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
   this->number_subproblems_solved++;
   return direction;
}

double SLP::compute_predicted_reduction(const Direction& direction, double step_length) const {
   // the predicted reduction is linear in the step length
   return -step_length * direction.objective;
}

int SLP::get_hessian_evaluation_count() const {
   // no second order evaluation is used
   return 0;
}