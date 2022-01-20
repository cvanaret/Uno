#include "QPSubproblem.hpp"
#include "solvers/QP/QPSolverFactory.hpp"

QPSubproblem::QPSubproblem(const Problem& problem, size_t max_number_variables, const Options& options) :
      Subproblem(problem.number_variables, max_number_variables, problem.number_constraints, false, NO_SOC, true,
            norm_from_string(options.at("residual_norm"))),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory::create(options.at("QP_solver"), max_number_variables, problem.number_constraints,
            problem.get_hessian_maximum_number_nonzeros() + max_number_variables, true)),
      // if no trust region is used, the problem should be convexified to guarantee boundedness + descent direction
      hessian_model(HessianModelFactory::create(options.at("hessian_model"), max_number_variables,
            problem.get_hessian_maximum_number_nonzeros() + max_number_variables, options.at("mechanism") != "TR", options)),
      initial_point(max_number_variables) {
}

void QPSubproblem::create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier,
      double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);

   // constraints
   current_iterate.evaluate_constraints(problem);

   // constraint Jacobian
   current_iterate.evaluate_constraint_jacobian(problem);
   this->constraint_jacobian = current_iterate.constraint_jacobian;

   // objective
   this->build_objective_model(problem, current_iterate, objective_multiplier);

   // bounds of the variables
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   // bounds of the linearized constraints
   this->set_constraints_bounds(problem, current_iterate.constraints);

   // reset the initial point
   initialize_vector(this->initial_point, 0.);
}

void QPSubproblem::build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // Hessian
   this->hessian_model->evaluate(problem, current_iterate.x, objective_multiplier, this->constraints_multipliers);

   // objective gradient
   this->set_scaled_objective_gradient(problem, current_iterate, objective_multiplier);

   // initial point
   initialize_vector(this->initial_point, 0.);
}

void QPSubproblem::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

Direction QPSubproblem::solve(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) {
   this->hessian_model->adjust_number_variables(this->number_variables);
   // compute QP direction
   Direction direction = this->solver->solve_QP(this->variables_bounds, this->constraints_bounds, this->objective_gradient,
         this->constraint_jacobian, *this->hessian_model->hessian, this->initial_point);
   this->number_subproblems_solved++;

   // compute dual *displacements* (note: SQP methods usually compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
   return direction;
}

PredictedReductionModel QPSubproblem::generate_predicted_reduction_model(const Problem& /*problem*/, const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() { // capture this and direction by reference
      // precompute expensive quantities
      const double linear_term = dot(direction.x, this->objective_gradient);
      const double quadratic_term = this->hessian_model->hessian->quadratic_product(direction.x, direction.x, this->number_variables) / 2.;
      // return a function of the step length that cheaply assembles the predicted reduction
      return [=](double step_length) { // capture the expensive quantities by value
         return -step_length * (linear_term + step_length * quadratic_term);
      };
   });
}

size_t QPSubproblem::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}