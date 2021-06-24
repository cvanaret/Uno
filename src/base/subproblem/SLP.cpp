#include "SLP.hpp"
#include "QPSolverFactory.hpp"

SLP::SLP(const Problem& problem, std::string QP_solver_name, bool /*use_trust_region*/) :
      Subproblem(problem),
      solver(QPSolverFactory::create(QP_solver_name, problem.number_variables, problem.number_constraints, 0, false)),
      initial_point(problem.number_variables) {
}

void SLP::generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);
   /* compute first- and second-order information */
   problem.evaluate_constraints(current_iterate.x, current_iterate.constraints);
   this->constraints_jacobian = problem.constraints_jacobian(current_iterate.x);

   this->objective_gradient = problem.objective_gradient(current_iterate.x);
   this->update_objective_multiplier(problem, current_iterate, objective_multiplier);

   /* bounds of the variables */
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->set_constraints_bounds(problem, current_iterate.constraints);

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

void SLP::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

Direction SLP::compute_direction(Statistics& /*statistics*/, const Problem& /*problem*/, Iterate& /*current_iterate*/) {
   /* solve the LP */
   Direction direction = this->solver->solve_LP(variables_bounds, constraints_bounds, this->objective_gradient,
         this->constraints_jacobian,this->initial_point);
   this->number_subproblems_solved++;
   DEBUG << direction;
   direction.predicted_reduction = [&](double step_length) {
      return SLP::compute_predicted_reduction_(direction, step_length);
   };
   return direction;
}

double SLP::compute_predicted_reduction_(Direction& direction, double step_length) {
   // the predicted reduction is linear in the step length
   return -step_length * direction.objective;
}

int SLP::get_hessian_evaluation_count() const {
   // no second order evaluation is used
   return 0;
}