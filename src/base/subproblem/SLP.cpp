#include "SLP.hpp"
#include "QPSolverFactory.hpp"

SLP::SLP(const Problem& problem, std::string QP_solver_name, bool /*use_trust_region*/, bool scale_residuals) :
      ActiveSetMethod(problem, scale_residuals),
      solver(QPSolverFactory::create(QP_solver_name, problem.number_variables, problem.number_constraints, 0, false)) {
}

void SLP::generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);
   /* compute first- and second-order information */
   problem.evaluate_constraints(current_iterate.x, this->constraints);
   this->constraints_jacobian = problem.constraints_jacobian(current_iterate.x);

   this->objective_gradient = problem.objective_gradient(current_iterate.x);
   this->update_objective_multiplier(problem, current_iterate, objective_multiplier);

   /* bounds of the variables */
   this->set_variables_bounds_(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->set_constraints_bounds(problem, this->constraints);

   /* set the initial point */
   clear(this->initial_point);
}

void SLP::update_objective_multiplier(const Problem& /*problem*/, const Iterate& /*current_iterate*/, double objective_multiplier) {
   // scale objective gradient
   if (objective_multiplier == 0.) {
      clear(this->objective_gradient);
   }
   else {
      if (objective_multiplier < 1.) {
         scale(this->objective_gradient, objective_multiplier);
      }
   }
}

Direction SLP::compute_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) {
   Direction direction = this->compute_lp_step_(problem, *this->solver, current_iterate, trust_region_radius);
   this->number_subproblems_solved++;
   DEBUG << direction;
   // attach the predicted reduction function
   direction.predicted_reduction = [&](double step_length) {
      return SLP::compute_predicted_reduction_(current_iterate, direction, step_length);
   };
   return direction;
}

double SLP::compute_predicted_reduction_(Iterate& current_iterate, Direction& direction, double step_length) {
   // the predicted reduction is quadratic in the step length
   if (step_length == 1.) {
      return -direction.objective;
   }
   else {
      return -step_length * dot(direction.x, current_iterate.objective_gradient);
   }
}

void SLP::evaluate_optimality_iterate_(const Problem& problem, Iterate& current_iterate) {
   /* compute first-order information */
   current_iterate.compute_objective_gradient(problem);
   current_iterate.compute_constraints_jacobian(problem);
}

void SLP::evaluate_feasibility_iterate_(const Problem& problem, Iterate& current_iterate, Direction& /*phase_2_direction*/) {
   /* compute first-order information */
   current_iterate.compute_constraints_jacobian(problem);
}

int SLP::get_hessian_evaluation_count() {
   return 0;
}