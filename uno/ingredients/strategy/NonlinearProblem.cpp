#include "NonlinearProblem.hpp"

NonlinearProblem::NonlinearProblem(const Model& model, size_t number_variables, size_t number_constraints):
      model(model), number_variables(number_variables), number_constraints(number_constraints) {
}

bool NonlinearProblem::is_constrained() const {
   return (0 < this->number_constraints);
}

void NonlinearProblem::evaluate_lagrangian_gradient(Iterate& iterate, std::vector<double>& lagrangian_gradient) const {
   iterate.evaluate_lagrangian_gradient(this->model, iterate.multipliers.constraints, iterate.multipliers.lower_bounds,
         iterate.multipliers.upper_bounds);
   lagrangian_gradient = iterate.lagrangian_gradient;
}

size_t NonlinearProblem::get_number_original_variables() const {
   return this->model.number_variables;
}

double NonlinearProblem::compute_constraint_lower_bound_violation(double constraint, size_t j) const {
   const double lower_bound = this->get_constraint_lower_bound(j);
   return std::max(0., lower_bound - constraint);
}

double NonlinearProblem::compute_constraint_upper_bound_violation(double constraint, size_t j) const {
   const double upper_bound = this->get_constraint_upper_bound(j);
   return std::max(0., constraint - upper_bound);
}

double NonlinearProblem::compute_constraint_violation(double constraint, size_t j) const {
   const double lower_bound_violation = this->compute_constraint_lower_bound_violation(constraint, j);
   const double upper_bound_violation = this->compute_constraint_upper_bound_violation(constraint, j);
   return std::max(lower_bound_violation, upper_bound_violation);
}

// compute ||c_S|| for a given set of constraints
double NonlinearProblem::compute_constraint_violation(const std::vector<double>& constraints, const std::vector<size_t>& constraint_set,
      Norm norm_type) const {
   auto residual_function = [&](size_t k) {
      const size_t j = constraint_set[k];
      return this->compute_constraint_violation(constraints[j], j);
   };
   return norm(residual_function, constraint_set.size(), norm_type);
}

// compute ||c||
double NonlinearProblem::compute_constraint_violation(const std::vector<double>& constraints, Norm norm_type) const {
   // create a lambda to avoid allocating an std::vector
   auto residual_function = [&](size_t j) {
      return this->compute_constraint_violation(constraints[j], j);
   };
   return norm(residual_function, constraints.size(), norm_type);
}


   /*
   const size_t number_original_variables = this->model.get_number_original_variables();
   initialize_vector(this->lagrangian_gradient, 0.);

   // objective gradient
   this->evaluate_objective_gradient(model);

   // scale the objective gradient with the objective multiplier
   this->original_evaluations.objective_gradient.for_each([&](size_t i, double derivative) {
      // in case there are additional variables, ignore them
      if (i < number_original_variables) {
         //this->lagrangian_gradient[i] += objective_multiplier * derivative;
         this->lagrangian_gradient[i] += derivative;
      }
   });

   // bound constraints
   for (size_t i = 0; i < number_original_variables; i++) {
      this->lagrangian_gradient[i] -= lower_bounds_multipliers[i] + upper_bounds_multipliers[i];
   }

   // constraints
   this->evaluate_constraint_jacobian(model);
   for (size_t j = 0; j < model.number_constraints; j++) {
      const double multiplier_j = constraint_multipliers[j];
      if (multiplier_j != 0.) {
         this->original_evaluations.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
            // in case there are additional variables, ignore them
            if (i < number_original_variables) {
               this->lagrangian_gradient[i] -= multiplier_j * derivative;
            }
         });
      }
   }
}
 */