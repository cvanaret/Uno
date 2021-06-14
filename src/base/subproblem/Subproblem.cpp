#include "Subproblem.hpp"
#include "LinearSolverFactory.hpp"

Subproblem::Subproblem(const Problem& problem, Norm residual_norm, bool scale_residuals) :
      number_variables(problem.number_variables),
      number_constraints(problem.number_constraints),
      variables_bounds(problem.number_variables),
      constraints_multipliers(problem.number_constraints),
      objective_value(INFINITY),
      // objective_gradient is a SparseVector
      constraints(problem.number_constraints),
      // constraints_jacobian is a vector of SparseVectors
      constraints_bounds(problem.number_constraints), initial_point(problem.number_variables), residual_norm(residual_norm),
      number_subproblems_solved(0), subproblem_definition_changed(false), scale_residuals(scale_residuals) {
}

/* compute least-square multipliers */
//if (0 < problem.number_constraints) {
//    first_iterate.compute_constraints_jacobian(problem);
//    first_iterate.multipliers.constraints = Subproblem::compute_least_square_multipliers(problem, first_iterate, multipliers.constraints, 1e4);
//}

Iterate Subproblem::evaluate_initial_point(const Problem& problem, const std::vector<double>& x, const Multipliers& multipliers) {
   Iterate first_iterate(x, multipliers);
   /* compute the optimality and feasibility measures of the initial point */
   this->compute_optimality_measures(problem, first_iterate);
   return first_iterate;
}

void Subproblem::compute_optimality_measures(const Problem& problem, Iterate& iterate) {
   iterate.compute_constraints(problem);
   // feasibility measure: residual of all constraints
   iterate.residuals.constraints = problem.compute_constraint_residual(iterate.constraints, this->residual_norm);
   // optimality
   iterate.compute_objective(problem);
   iterate.progress = {iterate.residuals.constraints, iterate.objective};
}


void Subproblem::project_point_in_bounds(std::vector<double>& x, const std::vector<Range>& variables_bounds) {
   for (size_t i = 0; i < x.size(); i++) {
      if (x[i] < variables_bounds[i].lb) {
         x[i] = variables_bounds[i].lb;
      }
      else if (variables_bounds[i].ub < x[i]) {
         x[i] = variables_bounds[i].ub;
      }
   }
}

double Subproblem::project_variable_in_interior(double variable_value, const Range& variable_bounds) {
   double k1 = 1e-2;
   double k2 = 1e-2;

   double perturbation_lb = std::min(k1 * std::max(1., std::abs(variable_bounds.lb)), k2 * (variable_bounds.ub - variable_bounds.lb));
   double perturbation_ub = std::min(k1 * std::max(1., std::abs(variable_bounds.ub)), k2 * (variable_bounds.ub - variable_bounds.lb));
   variable_value = std::max(variable_value, variable_bounds.lb + perturbation_lb);
   variable_value = std::min(variable_value, variable_bounds.ub - perturbation_ub);
   return variable_value;
}

void Subproblem::set_trust_region(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   this->set_variables_bounds_(problem, current_iterate, trust_region_radius);
}

void Subproblem::set_variables_bounds_(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   /* bounds intersected with trust region  */
   // very important: apply the trust region only on the original variables
   for (size_t i = 0; i < problem.number_variables; i++) {
      double lb = std::max(-trust_region_radius, problem.variables_bounds[i].lb - current_iterate.x[i]);
      double ub = std::min(trust_region_radius, problem.variables_bounds[i].ub - current_iterate.x[i]);
      this->variables_bounds[i] = {lb, ub};
   }
}

void Subproblem::set_constraints_bounds(const Problem& problem, const std::vector<double>& current_constraints) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      double lb = problem.constraint_bounds[j].lb - current_constraints[j];
      double ub = problem.constraint_bounds[j].ub - current_constraints[j];
      this->constraints_bounds[j] = {lb, ub};
   }
}

std::vector<double> Subproblem::compute_least_square_multipliers(const Problem& problem, Iterate& current_iterate,
      const std::vector<double>& default_multipliers, double multipliers_max_size) {
   std::unique_ptr<LinearSolver> linear_solver = LinearSolverFactory::create("MA57");
   return Subproblem::compute_least_square_multipliers(problem, current_iterate, default_multipliers, *linear_solver, multipliers_max_size);
}

std::vector<double> Subproblem::compute_least_square_multipliers(const Problem& problem, Iterate& current_iterate,
      const std::vector<double>& default_multipliers, LinearSolver& solver, double multipliers_max_size) {
   current_iterate.compute_objective_gradient(problem);
   current_iterate.compute_constraints_jacobian(problem);

   /******************************/
   /* build the symmetric matrix */
   /******************************/
   COOMatrix matrix(current_iterate.x.size() + problem.number_constraints, problem.hessian_maximum_number_nonzeros, 1);

   /* identity blocks */
   for (size_t i = 0; i < current_iterate.x.size(); i++) {
      matrix.insert(1., i, i);
   }
   /* Jacobian of general constraints */
   for (size_t j = 0; j < problem.number_constraints; j++) {
      for (const auto[variable_index, derivative]: current_iterate.constraints_jacobian[j]) {
         matrix.insert(derivative, variable_index, current_iterate.x.size() + j);
      }
   }
   DEBUG << "Multipliers estimation: KKT matrix:\n";
   for (size_t k = 0; k < matrix.matrix.size(); k++) {
      DEBUG << "m(" << matrix.row_indices[k] << ", " << matrix.column_indices[k] << ") = " << matrix.matrix[k] << "\n";
   }

   /********************************/
   /* generate the right-hand side */
   /********************************/
   std::vector<double> rhs(current_iterate.x.size() + problem.number_constraints);

   /* objective gradient */
   for (const auto[i, derivative]: current_iterate.objective_gradient) {
      rhs[i] += problem.objective_sign * derivative;
   }
   /* variable bound constraints */
   for (size_t i = 0; i < current_iterate.x.size(); i++) {
      rhs[i] -= current_iterate.multipliers.lower_bounds[i];
      rhs[i] -= current_iterate.multipliers.upper_bounds[i];
   }
   DEBUG << "Multipliers RHS:\n";
   print_vector(DEBUG, rhs);

   solver.do_symbolic_factorization(matrix);
   solver.solve(rhs);
   DEBUG << "Solution: ";
   std::vector<double>& solution = rhs;
   print_vector(DEBUG, solution);

   /* retrieve multipliers */
   std::vector<double> multipliers(problem.number_constraints);
   for (size_t j = 0; j < problem.number_constraints; j++) {
      multipliers[j] = solution[current_iterate.x.size() + j];
   }
   // if multipliers too big, discard them
   if (norm_inf(multipliers) > multipliers_max_size) {
      return default_multipliers;
   }
   return multipliers;
}

void Subproblem::compute_l1_linear_objective(const Iterate& current_iterate, const ConstraintPartition& constraint_partition) {
   /* objective function: sum of gradients of infeasible constraints */
   this->objective_gradient.clear();
   for (int j: constraint_partition.infeasible) {
      for (const auto [i, derivative]: current_iterate.constraints_jacobian[j]) {
         if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
            this->objective_gradient[i] -= derivative;
         }
         else {
            this->objective_gradient[i] += derivative;
         }
      }
   }
}

void Subproblem::generate_feasibility_bounds(const Problem& problem, const std::vector<double>& current_constraints, const ConstraintPartition&
constraint_partition) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      double lb, ub;
      if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
         lb = -INFINITY;
         ub = problem.constraint_bounds[j].lb - current_constraints[j];
      }
      else if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
         lb = problem.constraint_bounds[j].ub - current_constraints[j];
         ub = INFINITY;
      }
      else { // FEASIBLE
         lb = problem.constraint_bounds[j].lb - current_constraints[j];
         ub = problem.constraint_bounds[j].ub - current_constraints[j];
      }
      this->constraints_bounds[j] = {lb, ub};
   }
}

void Subproblem::compute_infeasibility_measures(const Problem& problem, Iterate& iterate, const ConstraintPartition& constraint_partition) {
   iterate.compute_constraints(problem);
   // feasibility measure: residual of all constraints
   double feasibility = problem.compute_constraint_residual(iterate.constraints, this->residual_norm);
   // optimality measure: residual of linearly infeasible constraints
   double objective = problem.compute_constraint_residual(iterate.constraints, constraint_partition.infeasible, this->residual_norm);
   iterate.progress = {feasibility, objective};
}

double Subproblem::compute_first_order_error(const Problem& problem, Iterate& iterate, double objective_multiplier) const {
   std::vector<double> lagrangian_gradient = iterate.lagrangian_gradient(problem, objective_multiplier, iterate.multipliers);
   return norm(lagrangian_gradient, this->residual_norm);
}

/* complementary slackness error. Use abs/1e-8 to safeguard */
double Subproblem::compute_complementarity_error_(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) const {
   double complementarity_error = 0.;
   /* bound constraints */
   for (size_t i = 0; i < iterate.x.size(); i++) {
      if (-INFINITY < problem.variables_bounds[i].lb) {
         complementarity_error += std::abs(multipliers.lower_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].lb));
      }
      if (problem.variables_bounds[i].ub < INFINITY) {
         complementarity_error += std::abs(multipliers.upper_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].ub));
      }
   }
   /* constraints */
   iterate.compute_constraints(problem);
   for (size_t j = 0; j < problem.number_constraints; j++) {
      double multiplier_j = multipliers.constraints[j];
      if (-INFINITY < problem.constraint_bounds[j].lb && 0. < multiplier_j) {
         complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].lb));
      }
      if (problem.constraint_bounds[j].ub < INFINITY && multiplier_j < 0.) {
         complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].ub));
      }
   }
   return complementarity_error;
}

void
Subproblem::compute_residuals(const Problem& problem, Iterate& iterate, const Multipliers& multipliers, double objective_multiplier) const {
   iterate.compute_constraints(problem);
   iterate.residuals.constraints = problem.compute_constraint_residual(iterate.constraints, this->residual_norm);
   // compute the KKT residual only if the objective multiplier is positive
   if (0. < objective_multiplier) {
      iterate.residuals.KKT = Subproblem::compute_first_order_error(problem, iterate, objective_multiplier);
   }
   else {
      iterate.residuals.KKT = Subproblem::compute_first_order_error(problem, iterate, 1.);
   }
   iterate.residuals.FJ = Subproblem::compute_first_order_error(problem, iterate, 0.);
   iterate.residuals.complementarity = this->compute_complementarity_error_(problem, iterate, multipliers);
   if (this->scale_residuals) {
      // TODO scale the residuals
   }
}
