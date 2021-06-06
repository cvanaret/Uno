#include <cmath>
#include "InteriorPoint.hpp"
#include "Uno.hpp"
#include "LinearSolverFactory.hpp"

InteriorPoint::InteriorPoint(Problem& problem, std::string linear_solver_name, std::string hessian_evaluation_method, bool use_trust_region,
      bool scale_residuals) : Subproblem(L2_NORM, problem.variables_bounds, scale_residuals), // use the l2 norm to compute residuals
/* if no trust region is used, the problem should be convexified. However, the inertia of the augmented matrix will be corrected later */
      hessian_evaluation(HessianEvaluationFactory::create(hessian_evaluation_method, problem.number_variables, problem
      .hessian_maximum_number_nonzeros, false)),
      linear_solver(LinearSolverFactory::create(linear_solver_name)), mu_optimality(0.1), mu_feasibility(mu_optimality),
      rhs_(problem.number_variables + problem.inequality_constraints.size() + problem.number_constraints), inertia_hessian_(0.),
      inertia_hessian_last_(0.), inertia_constraints_(0.), default_multiplier_(1.), iteration_(0), number_factorizations_(0),
      parameters_({0.99, 1e10, 100., 0.2, 1.5, 10., 1e10}) {
   /* the subproblem optimizes the original variables and slacks for inequality constraints */
   size_t number_variables = problem.number_variables + problem.inequality_constraints.size();
   this->bounds.resize(number_variables);

   /* identify the bounded variables */
   for (size_t i = 0; i < problem.number_variables; i++) {
      if (use_trust_region || (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
         this->lower_bounded_variables.insert(i);
      }
      if (use_trust_region || (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
         this->upper_bounded_variables.insert(i);
      }
   }
   /* identify the inequality constraint slacks */
   for (const auto[j, i]: problem.inequality_constraints) {
      int slack_index = problem.number_variables + i;
      if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
         this->lower_bounded_variables.insert(slack_index);
      }
      if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
         this->upper_bounded_variables.insert(slack_index);
      }
      // register the bounds of the slacks
      this->bounds[slack_index] = problem.constraint_bounds[j];
   }
}

Iterate
InteriorPoint::evaluate_initial_point(const Problem& problem, const std::vector<double>& x, const Multipliers& default_multipliers) {
   int number_variables = problem.number_variables + problem.inequality_constraints.size();

   /* make the initial point strictly feasible */
   std::vector<double> reformulated_x(number_variables);
   for (size_t i = 0; i < problem.number_variables; i++) {
      reformulated_x[i] = Subproblem::project_strictly_variable_in_bounds(x[i], problem.variables_bounds[i]);
   }

   Multipliers multipliers(number_variables, problem.number_constraints);
   /* generate the bound multipliers */
   for (int i: this->lower_bounded_variables) {
      multipliers.lower_bounds[i] = this->default_multiplier_; // positive multiplier
   }
   for (int i: this->upper_bounded_variables) {
      multipliers.upper_bounds[i] = -this->default_multiplier_; // negative multiplier
   }

   /* generate the first iterate */
   Iterate first_iterate(reformulated_x, multipliers);

   /* initialize the slacks */
   first_iterate.compute_constraints(problem);
   for (const auto[j, i]: problem.inequality_constraints) {
      int slack_index = problem.number_variables + i;
      double slack_value = Subproblem::project_strictly_variable_in_bounds(first_iterate.constraints[j], problem.constraint_bounds[j]);
      first_iterate.x[slack_index] = slack_value;
   }

   /* evaluate the constraint Jacobian */
   first_iterate.compute_constraints_jacobian(problem);
   // contribution of the slacks
   for (const auto[j, i]: problem.inequality_constraints) {
      int slack_index = problem.number_variables + i;
      first_iterate.constraints_jacobian[j][slack_index] = -1.;
   }
   /* compute least-square multipliers */
   if (0 < problem.number_constraints) {
      first_iterate.multipliers.constraints =
            Subproblem::compute_least_square_multipliers(problem, first_iterate, default_multipliers.constraints, *this->linear_solver);
   }

   DEBUG << problem.inequality_constraints.size() << " slacks\n";
   DEBUG << first_iterate.multipliers.lower_bounds.size() << " bound multipliers\n";
   DEBUG << first_iterate.multipliers.constraints.size() << " constraint multipliers\n";
   DEBUG << "variable lb: ";
   print_vector(DEBUG, this->lower_bounded_variables);
   DEBUG << "variable ub: ";
   print_vector(DEBUG, this->upper_bounded_variables);

   /* compute the optimality and feasibility measures of the initial point */
   this->compute_optimality_measures(problem, first_iterate);
   return first_iterate;
}

double InteriorPoint::compute_KKT_error_scaling_(Iterate& current_iterate) {
   /* KKT error */
   double norm_1_constraint_multipliers = norm_1(current_iterate.multipliers.constraints);
   double norm_1_bound_multipliers = norm_1(current_iterate.multipliers.lower_bounds) + norm_1(current_iterate.multipliers.upper_bounds);
   double sd = std::max(this->parameters_.smax, (norm_1_constraint_multipliers + norm_1_bound_multipliers) /
                                                (current_iterate.x.size() + current_iterate.multipliers.constraints.size())) /
               this->parameters_.smax;
   return sd;
}

/* reduced primal-dual approach */
std::vector<Direction>
InteriorPoint::compute_directions(Problem& problem, Iterate& current_iterate, double /*objective_multiplier*/, double trust_region_radius) {
   DEBUG << "\nCurrent iterate: " << current_iterate;

   current_iterate.compute_constraints_jacobian(problem);
   // contribution of the slacks
   for (const auto[j, i]: problem.inequality_constraints) {
      int slack_index = problem.number_variables + i;
      current_iterate.constraints_jacobian[j][slack_index] = -1.;
   }

   /* scaled error terms */
   this->compute_residuals(problem, current_iterate, current_iterate.multipliers, 1.);
   double sd = this->compute_KKT_error_scaling_(current_iterate);
   double KKTerror = current_iterate.residuals.KKT / sd;
   double central_complementarity_error = this->compute_central_complementarity_error(current_iterate, this->mu_optimality, this->bounds);
   DEBUG << "IPM error (KKT: " << KKTerror << ", cmpl: " << central_complementarity_error << ", feas: "
         << current_iterate.residuals.constraints << ")\n";

   /* update of the barrier problem */
   double error = std::max(KKTerror, std::max(central_complementarity_error, current_iterate.residuals.constraints));
   if (error <= this->parameters_.k_epsilon * this->mu_optimality) {
      // TODO pass tolerance
      double tolerance = 1e-8;
      this->mu_optimality = std::max(tolerance / 10.,
            std::min(this->parameters_.k_mu * this->mu_optimality, std::pow(this->mu_optimality, this->parameters_.theta_mu)));
      DEBUG << "IPM: mu updated to " << this->mu_optimality << " and filter reset\n";
      // signal the redefinition of the problem to the globalization strategy
      this->subproblem_definition_changed = true;
   }
   DEBUG << "mu is " << this->mu_optimality << "\n";
   this->iteration_++;

   /* evaluate the functions at the current iterate */
   this->evaluate_optimality_iterate_(problem, current_iterate);

   /************************/
   /* solve the KKT system */
   /************************/
   try {
      /* KKT matrix */
      COOMatrix kkt_matrix = this->generate_optimality_kkt_matrix_(problem, current_iterate, this->bounds);

      /* inertia correction (includes factorization) */
      this->modify_inertia_(kkt_matrix, current_iterate.x.size(), problem.number_constraints, problem.type);
      DEBUG << "KKT matrix:\n" << kkt_matrix << "\n";

      /* right-hand side */
      this->generate_kkt_rhs_(problem, current_iterate);

      /* compute the solution (Δx, -Δλ) */
      this->linear_solver->solve(this->rhs_);
      this->number_subproblems_solved++;
      std::vector<double>& solution_IPM = this->rhs_;

      /* generate IPM direction */
      Direction direction = this->generate_direction_(problem, current_iterate, solution_IPM);
      direction.status = OPTIMAL;
      direction.phase = OPTIMALITY;
      direction.norm = norm_inf(direction.x, problem.number_variables);
      direction.predicted_reduction = this->compute_predicted_reduction_;
      /* evaluate the barrier objective */
      direction.objective = this->evaluate_local_model_(problem, current_iterate, direction.x);
      return std::vector<Direction>{direction};
   }
   catch (const UnstableInertiaCorrection& e) {
      /* unstable factorization during optimality phase */
      Direction direction(current_iterate.x, current_iterate.multipliers); // dummy solution
      return this->restore_feasibility(problem, current_iterate, direction, trust_region_radius);
   }
}

void InteriorPoint::evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate) {
   /* compute barrier gradient */
   current_iterate.compute_objective_gradient(problem);
   // contribution of bound constraints
   for (int i: this->lower_bounded_variables) {
      current_iterate.objective_gradient[i] -= this->mu_optimality / (current_iterate.x[i] - this->bounds[i].lb);
   }
   for (int i: this->upper_bounded_variables) {
      current_iterate.objective_gradient[i] -= this->mu_optimality / (current_iterate.x[i] - this->bounds[i].ub);
   }

   /* compute constraint Jacobian */
   current_iterate.compute_constraints_jacobian(problem);
   // contribution of the slacks
   for (const auto[j, i]: problem.inequality_constraints) {
      int slack_index = problem.number_variables + i;
      current_iterate.constraints_jacobian[j][slack_index] = -1.;
   }

   /* compute second-order information */
   this->hessian_evaluation->compute(problem, current_iterate.x, problem.objective_sign, current_iterate.multipliers.constraints);
}

Direction InteriorPoint::generate_direction_(Problem& problem, Iterate& current_iterate, std::vector<double>& solution_IPM) {
   int number_variables = problem.number_variables + problem.inequality_constraints.size();

   /* retrieve +Δλ (Nocedal p590) */
   for (size_t j = 0; j < problem.number_constraints; j++) {
      int multiplier_index = number_variables + j;
      solution_IPM[multiplier_index] = -solution_IPM[multiplier_index];
   }

   /* compute bound multiplier displacements Δz */
   std::vector<double>
         lower_delta_z = this->compute_lower_bound_multiplier_displacements_(current_iterate, solution_IPM, bounds, this->mu_optimality);
   std::vector<double>
         upper_delta_z = this->compute_upper_bound_multiplier_displacements_(current_iterate, solution_IPM, bounds, this->mu_optimality);

   /* "fraction to boundary" rule for variables and bound multipliers */
   std::vector<double> trial_x(current_iterate.x.size());
   Multipliers trial_multipliers(current_iterate.x.size(), problem.number_constraints);
   double tau = std::max(this->parameters_.tau_min, 1. - this->mu_optimality);
   // scale primal variables and constraints multipliers
   double primal_length = this->compute_primal_length_(current_iterate, solution_IPM, bounds, tau);
   for (int i = 0; i < number_variables; i++) {
      trial_x[i] = primal_length * solution_IPM[i];
   }
   for (size_t j = 0; j < problem.number_constraints; j++) {
      trial_multipliers.constraints[j] = current_iterate.multipliers.constraints[j] + primal_length * solution_IPM[number_variables + j];
   }

   // scale dual variables
   double dual_length = this->compute_dual_length_(current_iterate, tau, lower_delta_z, upper_delta_z);
   for (int i = 0; i < number_variables; i++) {
      trial_multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_length * lower_delta_z[i];
      trial_multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_length * upper_delta_z[i];
      // rescale the multipliers (IPOPT paper p6)
      //trial_multipliers.lower_bounds[i] = std::max(std::min(trial_multipliers.lower_bounds[i], this->parameters.kappa*this->mu_optimality/(trial_x[i] - subproblem_variables_bounds[i].lb)), this->mu_optimality/(this->parameters.kappa*(trial_x[i] - subproblem_variables_bounds[i].lb)));
   }

   DEBUG << "MA57 solution:\n";
   DEBUG << "Δx: ";
   print_vector(DEBUG, solution_IPM, 0, problem.number_variables);
   DEBUG << "Δs: ";
   print_vector(DEBUG, solution_IPM, problem.number_variables, problem.inequality_constraints.size());
   DEBUG << "Δz_L: ";
   print_vector(DEBUG, lower_delta_z);
   DEBUG << "Δz_U: ";
   print_vector(DEBUG, upper_delta_z);
   DEBUG << "Δλ: ";
   print_vector(DEBUG, solution_IPM, number_variables, problem.number_constraints);
   DEBUG << "primal length = " << primal_length << "\n";
   DEBUG << "dual length = " << dual_length << "\n\n";

   Direction direction(trial_x, trial_multipliers);
   return direction;
}

double
InteriorPoint::compute_primal_length_(Iterate& current_iterate, std::vector<double>& ipm_solution, std::vector<Range>& variables_bounds,
      double tau) {
   double primal_length = 1.;
   for (int i: this->lower_bounded_variables) {
      double trial_alpha_xi = -tau * (current_iterate.x[i] - variables_bounds[i].lb) / ipm_solution[i];
      if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   for (int i: this->upper_bounded_variables) {
      double trial_alpha_xi = -tau * (current_iterate.x[i] - variables_bounds[i].ub) / ipm_solution[i];
      if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   return primal_length;
}

double InteriorPoint::compute_dual_length_(Iterate& current_iterate, double tau, std::vector<double>& lower_delta_z,
      std::vector<double>& upper_delta_z) {
   double dual_length = 1.;
   for (size_t i = 0; i < current_iterate.multipliers.lower_bounds.size(); i++) {
      double trial_alpha_zj = -tau * current_iterate.multipliers.lower_bounds[i] / lower_delta_z[i];
      if (0 < trial_alpha_zj && trial_alpha_zj <= 1.) {
         dual_length = std::min(dual_length, trial_alpha_zj);
      }
      trial_alpha_zj = -tau * current_iterate.multipliers.upper_bounds[i] / upper_delta_z[i];
      if (0 < trial_alpha_zj && trial_alpha_zj <= 1.) {
         dual_length = std::min(dual_length, trial_alpha_zj);
      }
   }
   return dual_length;
}

COOMatrix InteriorPoint::generate_optimality_kkt_matrix_(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
   int number_variables = problem.number_variables + problem.inequality_constraints.size();

   /* compute the Lagrangian Hessian */
   COOMatrix kkt_matrix = this->hessian_evaluation->hessian.to_COO();
   kkt_matrix.dimension = number_variables + problem.number_constraints;

   /* variable bound constraints */
   for (int i: this->lower_bounded_variables) {
      kkt_matrix.insert(current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - variables_bounds[i].lb), i, i);
   }
   for (int i: this->upper_bounded_variables) {
      kkt_matrix.insert(current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - variables_bounds[i].ub), i, i);
   }

   /* Jacobian of general constraints */
   for (size_t j = 0; j < problem.number_constraints; j++) {
      for (const auto[variable_index, derivative]: current_iterate.constraints_jacobian[j]) {
         kkt_matrix.insert(derivative, variable_index, number_variables + j);
      }
   }
   return kkt_matrix;
}

void InteriorPoint::factorize_(COOMatrix& kkt_matrix, FunctionType problem_type) {
   // compute the symbolic factorization only when:
   // the problem has a non constant Hessian (ie is not an LP or a QP) or it is the first factorization
   // TODO: for QPs as well, but only when the sparsity pattern is constant
   if (force_symbolic_factorization || problem_type == LINEAR || this->number_factorizations_ == 0) {
      this->linear_solver->do_symbolic_factorization(kkt_matrix);
   }
   this->linear_solver->do_numerical_factorization(kkt_matrix);
   this->number_factorizations_++;
}

void InteriorPoint::modify_inertia_(COOMatrix& kkt_matrix, int size_first_block, int size_second_block, FunctionType problem_type) {
   this->inertia_hessian_ = 0.;
   this->inertia_constraints_ = 0.;
   DEBUG << "Testing factorization with inertia term " << this->inertia_hessian_ << "\n";
   this->factorize_(kkt_matrix, problem_type);

   bool good_inertia = false;
   if (!this->linear_solver->matrix_is_singular() && this->linear_solver->number_negative_eigenvalues() == size_second_block) {
      DEBUG << "Factorization was a success\n";
      good_inertia = true;
   }
   else {
      // inertia term for constraints
      if (this->linear_solver->matrix_is_singular()) {
         DEBUG << "Matrix is singular\n";
         this->inertia_constraints_ = 1e-8 * std::pow(this->mu_optimality, 0.25);
      }
      else {
         this->inertia_constraints_ = 0.;
      }
      // inertia term for Hessian
      if (this->inertia_hessian_last_ == 0.) {
         this->inertia_hessian_ = 1e-4;
      }
      else {
         this->inertia_hessian_ = std::max(1e-20, this->inertia_hessian_last_ / 3.);
      }
   }

   int current_matrix_size = kkt_matrix.matrix.size();
   if (!good_inertia) {
      for (int i = 0; i < size_first_block; i++) {
         kkt_matrix.insert(this->inertia_hessian_, i, i);
      }
      for (int j = size_first_block; j < size_first_block + size_second_block; j++) {
         kkt_matrix.insert(-this->inertia_constraints_, j, j);
      }
   }

   while (!good_inertia) {
      DEBUG << "Testing factorization with inertia term " << this->inertia_hessian_ << "\n";
      this->factorize_(kkt_matrix, problem_type);

      if (!this->linear_solver->matrix_is_singular() && this->linear_solver->number_negative_eigenvalues() == size_second_block) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n";
         this->inertia_hessian_last_ = this->inertia_hessian_;
      }
      else {
         if (this->inertia_hessian_last_ == 0.) {
            this->inertia_hessian_ *= 100.;
         }
         else {
            this->inertia_hessian_ *= 8.;
         }
         if (1e40 < this->inertia_hessian_) {
            throw UnstableInertiaCorrection();
         }
         else {
            for (int i = 0; i < size_first_block; i++) {
               kkt_matrix.matrix[current_matrix_size + i] = this->inertia_hessian_;
            }
            for (int j = size_first_block; j < size_first_block + size_second_block; j++) {
               kkt_matrix.matrix[current_matrix_size + j] = -this->inertia_constraints_;
            }
         }
      }
   }
}

void InteriorPoint::generate_kkt_rhs_(Problem& problem, Iterate& current_iterate) {
   int number_variables = problem.number_variables + problem.inequality_constraints.size();

   /* generate the right-hand side */
   for (size_t i = 0; i < this->rhs_.size(); i++) {
      this->rhs_[i] = 0.;
   }

   /* barrier objective gradient */
   for (const auto[i, derivative]: current_iterate.objective_gradient) {
      this->rhs_[i] = -derivative;
   }

   /* constraint gradients */
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (current_iterate.multipliers.constraints[j] != 0.) {
         for (const auto[variable_index, derivative]: current_iterate.constraints_jacobian[j]) {
            this->rhs_[variable_index] += current_iterate.multipliers.constraints[j] * derivative;
         }
      }
   }

   /* constraint evaluations */
   int slack_index = 0;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (problem.constraint_status[j] == EQUAL_BOUNDS) {
         // add the bound
         this->rhs_[number_variables + j] = -(current_iterate.constraints[j] - problem.constraint_bounds[j].lb);
      }
      else {
         // add the slack
         this->rhs_[number_variables + j] = -(current_iterate.constraints[j] - current_iterate.x[problem.number_variables + slack_index]);
         slack_index++;
      }
   }
   DEBUG << "RHS: ";
   print_vector(DEBUG, this->rhs_);
}

std::vector<double> InteriorPoint::compute_lower_bound_multiplier_displacements_(Iterate& current_iterate, std::vector<double>& solution,
      std::vector<Range>& variables_bounds, double mu) {
   std::vector<double> delta_z(current_iterate.multipliers.lower_bounds.size());
   for (int i: this->lower_bounded_variables) {
      delta_z[i] = mu / (current_iterate.x[i] - variables_bounds[i].lb) - current_iterate.multipliers.lower_bounds[i] -
                   current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - variables_bounds[i].lb) * solution[i];
   }
   return delta_z;
}

std::vector<double> InteriorPoint::compute_upper_bound_multiplier_displacements_(Iterate& current_iterate, std::vector<double>& solution,
      std::vector<Range>& variables_bounds, double mu) {
   std::vector<double> delta_z(current_iterate.multipliers.upper_bounds.size());
   for (int i: this->upper_bounded_variables) {
      delta_z[i] = mu / (current_iterate.x[i] - variables_bounds[i].ub) - current_iterate.multipliers.upper_bounds[i] -
                   current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - variables_bounds[i].ub) * solution[i];
   }
   return delta_z;
}

void InteriorPoint::compute_optimality_measures(const Problem& problem, Iterate& iterate) {
   /* evaluate constraints with slacks */
   iterate.feasibility_measure = this->constraint_violation(problem, iterate);
   /* compute barrier objective */
   iterate.optimality_measure = this->barrier_function_(problem, iterate, this->bounds);
}

void InteriorPoint::compute_infeasibility_measures(const Problem& problem, Iterate& iterate, const Direction& /*direction*/) {
   this->compute_optimality_measures(problem, iterate);
}

double InteriorPoint::constraint_violation(const Problem& problem, Iterate& iterate) {
   iterate.compute_constraints(problem);
   // compute l2 square norm
   std::vector<double> residuals(problem.number_constraints);
   size_t slack_index = problem.number_variables;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (problem.constraint_status[j] == EQUAL_BOUNDS) {
         double constraint_value = iterate.constraints[j] - problem.constraint_bounds[j].lb;
         residuals[j] = constraint_value;
      }
      else {
         double constraint_value = iterate.constraints[j] - iterate.x[slack_index];
         residuals[j] = constraint_value;
         slack_index++;
      }
   }
   return norm(residuals, this->residual_norm);
}

double InteriorPoint::barrier_function_(const Problem& problem, Iterate& iterate, const std::vector<Range>& variables_bounds) {
   /* original objective */
   iterate.compute_objective(problem);
   double objective = iterate.objective;

   /* bound constraints */
   for (int i: this->lower_bounded_variables) {
      objective -= this->mu_optimality * std::log(iterate.x[i] - variables_bounds[i].lb);
   }
   for (int i: this->upper_bounded_variables) {
      objective -= this->mu_optimality * std::log(variables_bounds[i].ub - iterate.x[i]);
   }
   return objective;
}

double InteriorPoint::evaluate_local_model_(Problem& /*problem*/, Iterate& current_iterate, std::vector<double>& solution) {
   double subproblem_objective = dot(solution, current_iterate.objective_gradient);
   return subproblem_objective;
}

double
InteriorPoint::compute_predicted_reduction_(Problem& /*problem*/, Iterate& /*current_iterate*/, Direction& direction, double step_length) {
   // the predicted reduction is linear
   return -step_length * direction.objective;
}

std::vector<Direction> InteriorPoint::restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& /*phase_2_direction*/,
      double /*trust_region_radius*/) {
   int number_variables = problem.number_variables + problem.inequality_constraints.size();

   DEBUG << "restoration x: ";
   print_vector(DEBUG, current_iterate.x);

   /* multipliers = 2*c */
   std::vector<double> restoration_multipliers(problem.number_constraints);
   int slack_index = 0;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (problem.constraint_status[j] == EQUAL_BOUNDS) {
         // add the bound
         restoration_multipliers[j] = 2 * (current_iterate.constraints[j] - problem.constraint_bounds[j].lb);
      }
      else {
         // add the slack
         restoration_multipliers[j] = 2 * (current_iterate.constraints[j] - current_iterate.x[problem.number_variables + slack_index]);
         slack_index++;
      }
   }

   /* compute the Lagrangian Hessian */
   this->hessian_evaluation->compute(problem, current_iterate.x, 0., restoration_multipliers);
   UnoMatrix kkt_matrix = this->hessian_evaluation->hessian.to_UnoMatrix(number_variables);
   // contribution of 2 \nabla c \nabla c^T
   for (size_t j = 0; j < problem.number_constraints; j++) {
      DEBUG << "Gradient c" << j << ": ";
      print_vector(DEBUG, current_iterate.constraints_jacobian[j]);
      kkt_matrix.add_outer_product(current_iterate.constraints_jacobian[j], 2.);
   }
   // variable bound constraints
   for (int i: this->lower_bounded_variables) {
      kkt_matrix.insert(current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - this->bounds[i].lb), i, i);
   }
   for (int i: this->upper_bounded_variables) {
      kkt_matrix.insert(current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - this->bounds[i].ub), i, i);
   }

   /* factorization by the linear solver */
   COOMatrix coo_matrix = kkt_matrix.to_COO();

   /* inertia correction */
   this->modify_inertia_(coo_matrix, current_iterate.x.size(), 0, problem.type);

   DEBUG << "restoration KKT matrix:\n" << coo_matrix;

   /* right-hand side */
   std::vector<double> rhs(number_variables);
   // constraint Jacobian
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (restoration_multipliers[j] != 0.) {
         for (const auto[i, derivative]: current_iterate.constraints_jacobian[j]) {
            rhs[i] += restoration_multipliers[j] * derivative;
         }
      }
   }
   // variable bound constraints
   for (int i: this->lower_bounded_variables) {
      rhs[i] += this->mu_feasibility / (current_iterate.x[i] - this->bounds[i].lb);
   }
   for (int i: this->upper_bounded_variables) {
      rhs[i] += this->mu_feasibility / (current_iterate.x[i] - this->bounds[i].ub);
   }
   DEBUG << "restoration RHS: ";
   print_vector(DEBUG, rhs);
   DEBUG << "\n";

   /* compute the solution Δx */
   this->linear_solver->solve(rhs);
   this->number_subproblems_solved++;
   std::vector<double>& solution_IPM = rhs;

   /* compute bound multiplier displacements Δz */
   std::vector<double>
         lower_delta_z = this->compute_lower_bound_multiplier_displacements_(current_iterate, solution_IPM, bounds, this->mu_feasibility);
   std::vector<double>
         upper_delta_z = this->compute_upper_bound_multiplier_displacements_(current_iterate, solution_IPM, bounds, this->mu_feasibility);

   /* create the solution */
   std::vector<double> trial_x(current_iterate.x.size());
   Multipliers trial_multipliers(current_iterate.x.size(), current_iterate.constraints.size());
   double tau = std::max(this->parameters_.tau_min, 1. - this->mu_feasibility);
   // scale primal variables and constraints multipliers
   double primal_length = this->compute_primal_length_(current_iterate, solution_IPM, bounds, tau);
   for (int i = 0; i < number_variables; i++) {
      trial_x[i] = primal_length * solution_IPM[i];
   }
   // scale dual variables
   double dual_length = this->compute_dual_length_(current_iterate, tau, lower_delta_z, upper_delta_z);
   for (size_t i = 0; i < current_iterate.multipliers.lower_bounds.size(); i++) {
      trial_multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_length * lower_delta_z[i];
      trial_multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_length * upper_delta_z[i];
      // TODO rescale the multipliers (IPOPT paper p6)
   }

   DEBUG << "MA57 restoration solution:\n";
   DEBUG << "Δx: ";
   print_vector(DEBUG, solution_IPM, 0, problem.number_variables);
   DEBUG << "Δs: ";
   print_vector(DEBUG, solution_IPM, problem.number_variables, problem.inequality_constraints.size());
   DEBUG << "Δz_L: ";
   print_vector(DEBUG, lower_delta_z);
   DEBUG << "Δz_U: ";
   print_vector(DEBUG, upper_delta_z);
   DEBUG << "primal length = " << primal_length << "\n";
   DEBUG << "dual length = " << dual_length << "\n\n";

   Direction direction(trial_x, trial_multipliers);
   direction.status = INFEASIBLE;
   direction.phase = RESTORATION;
   direction.norm = norm_inf(direction.x, problem.number_variables);
//    direction.predicted_reduction = [&](double step_length) {
//        return this->compute_predicted_reduction_(direction, step_length);
//    };
   direction.predicted_reduction = this->compute_predicted_reduction_;
   return std::vector<Direction>{direction};
}

double InteriorPoint::compute_central_complementarity_error(Iterate& iterate, double mu, std::vector<Range>& variables_bounds) {
   std::vector<double> residuals(iterate.x.size());
   /* variable bound constraints */
   for (size_t i = 0; i < iterate.x.size(); i++) {
      if (-INFINITY < variables_bounds[i].lb) {
         residuals[i] = iterate.multipliers.lower_bounds[i] * (iterate.x[i] - variables_bounds[i].lb) - mu;
      }
      if (variables_bounds[i].ub < INFINITY) {
         residuals[i] = iterate.multipliers.upper_bounds[i] * (iterate.x[i] - variables_bounds[i].ub) - mu;
      }
   }

   /* scaling */
   double sc = std::max(this->parameters_.smax,
         (norm_1(iterate.multipliers.lower_bounds) + norm_1(iterate.multipliers.upper_bounds)) / iterate.x.size()) / this->parameters_.smax;
   return norm(residuals, this->residual_norm) / sc;
}
