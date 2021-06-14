#include "SQP.hpp"
#include "QPSolverFactory.hpp"

SQP::SQP(const Problem& problem, const std::string& QP_solver_name, const std::string& hessian_evaluation_method, bool use_trust_region,
      bool scale_residuals) : ActiveSetMethod(problem, scale_residuals),
      solver(QPSolverFactory::create(QP_solver_name, problem.number_variables, problem.number_constraints,
            // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
            problem.hessian_maximum_number_nonzeros + problem.number_variables, true)),
      /* if no trust region is used, the problem should be convexified by controlling the inertia of the Hessian */
      hessian_evaluation(HessianEvaluationFactory::create(hessian_evaluation_method, problem.number_variables, problem.hessian_maximum_number_nonzeros,
          !use_trust_region)) {
}

void SQP::generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
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

void SQP::update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) {
   // evaluate the Hessian
   this->hessian_evaluation->compute(problem, current_iterate.x, objective_multiplier, this->constraints_multipliers);

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

Direction SQP::compute_direction(const Problem& /*problem*/, Iterate& current_iterate, double /*trust_region_radius*/) {
   /* compute QP direction */
   Direction direction = this->solver->solve_QP(this->variables_bounds, constraints_bounds, this->objective_gradient,
         this->constraints_jacobian, this->hessian_evaluation->hessian, this->initial_point);
   this->number_subproblems_solved++;
   DEBUG << direction;
   // attach the predicted reduction function
   direction.predicted_reduction = [&](double step_length) {
      return this->compute_predicted_reduction_(current_iterate, direction, step_length);
   };
   return direction;
}

double SQP::compute_predicted_reduction_(Iterate& current_iterate, Direction& direction, double step_length) const {
   // the predicted reduction is quadratic in the step length
   if (step_length == 1.) {
      return -direction.objective;
   }
   else {
      double linear_term = dot(direction.x, current_iterate.objective_gradient);
      double quadratic_term = this->hessian_evaluation->hessian.quadratic_product(direction.x, direction.x) / 2.;
      return -step_length * (linear_term + step_length * quadratic_term);
   }
}

int SQP::get_hessian_evaluation_count() {
   return this->hessian_evaluation->evaluation_count;
}