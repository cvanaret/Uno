#include <cmath>
#include "InteriorPoint.hpp"
#include "LinearSolverFactory.hpp"

InteriorPoint::InteriorPoint(const Problem& problem, const std::string& linear_solver_name, const std::string& hessian_evaluation_method, bool
use_trust_region, bool scale_residuals) :
      Subproblem(problem, scale_residuals), // use the l2 norm to compute residuals
      barrier_parameter(0.1),
/* if no trust region is used, the problem should be convexified. However, the inertia of the augmented matrix will be corrected later */
      hessian_evaluation(HessianEvaluationFactory::create(hessian_evaluation_method, problem.number_variables, problem
      .hessian_maximum_number_nonzeros, false)),
      linear_solver(LinearSolverFactory::create(linear_solver_name)),
      parameters({0.99, 1e10, 100., 0.2, 1.5, 10., 1e10}),
      rhs(problem.number_variables + problem.inequality_constraints.size() + problem.number_constraints), inertia_hessian(0.),
      inertia_hessian_last_(0.), inertia_constraints(0.), default_multiplier_(1.), iteration(0), number_factorizations_(0) {

   /* the subproblem optimizes the original variables and slacks for inequality constraints */
   this->number_variables = problem.number_variables + problem.inequality_constraints.size();
   this->variables_bounds.resize(this->number_variables);

   for (size_t i = 0; i < problem.number_variables; i++) {
      this->variables_bounds[i] = problem.variables_bounds[i];
   }

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
      size_t slack_index = problem.number_variables + i;
      if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
         this->lower_bounded_variables.insert(slack_index);
      }
      if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
         this->upper_bounded_variables.insert(slack_index);
      }
      // register the bounds of the slacks
      this->variables_bounds[slack_index] = problem.constraint_bounds[j];
   }
}

Iterate InteriorPoint::generate_initial_point(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   statistics.add_column("barrier param.", Statistics::double_width, 8);

   // resize to the new size (primals + slacks)
   x.resize(this->number_variables);
   multipliers.lower_bounds.resize(this->number_variables);
   multipliers.upper_bounds.resize(this->number_variables);

   /* make the initial point strictly feasible */
   for (size_t i = 0; i < problem.number_variables; i++) {
      x[i] = Subproblem::push_variable_to_interior(x[i], problem.variables_bounds[i]);
   }

   /* set the bound multipliers */
   for (size_t i: this->lower_bounded_variables) {
      multipliers.lower_bounds[i] = this->default_multiplier_; // positive multiplier
   }
   for (size_t i: this->upper_bounded_variables) {
      multipliers.upper_bounds[i] = -this->default_multiplier_; // negative multiplier
   }

   /* generate the first iterate */
   Iterate first_iterate(x, multipliers);

   /* initialize the slacks and add contribution to the constraint Jacobian */
   first_iterate.compute_constraints(problem);
   first_iterate.compute_constraints_jacobian(problem);
   for (const auto[j, i]: problem.inequality_constraints) {
      double slack_value = Subproblem::push_variable_to_interior(first_iterate.constraints[j], problem.constraint_bounds[j]);
      first_iterate.x[problem.number_variables + i] = slack_value;
      first_iterate.constraints_jacobian[j][problem.number_variables + i] = -1.;
   }

   /* compute least-square multipliers */
   if (0 < problem.number_constraints) {
      Subproblem::compute_least_square_multipliers(problem, first_iterate, first_iterate.multipliers.constraints);
   }

   DEBUG << problem.inequality_constraints.size() << " slacks\n";
   DEBUG << first_iterate.multipliers.lower_bounds.size() << " bound multipliers\n";
   DEBUG << first_iterate.multipliers.constraints.size() << " constraint multipliers\n";
   DEBUG << "variable lb: ";
   print_vector(DEBUG, this->lower_bounded_variables);
   DEBUG << "variable ub: ";
   print_vector(DEBUG, this->upper_bounded_variables);

   /* compute the optimality and feasibility measures of the initial point */
   this->compute_progress_measures(problem, first_iterate);
   return first_iterate;
}

void InteriorPoint::set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   /* here, we work with the original bounds  */
   // very important: apply the trust region only on the original variables (not the slacks)
   for (size_t i = 0; i < problem.number_variables; i++) {
      double lb = std::max(current_iterate.x[i] - trust_region_radius, problem.variables_bounds[i].lb);
      double ub = std::min(current_iterate.x[i] + trust_region_radius, problem.variables_bounds[i].ub);
      this->variables_bounds[i] = {lb, ub};
   }
}

void InteriorPoint::generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);
   /* compute first- and second-order information */
   problem.evaluate_constraints(current_iterate.x, this->constraints);
   this->constraints_jacobian = problem.constraints_jacobian(current_iterate.x);
   // contribution of the slacks
   for (const auto[j, i]: problem.inequality_constraints) {
      this->constraints_jacobian[j][problem.number_variables + i] = -1.;
   }

   this->objective_gradient = problem.objective_gradient(current_iterate.x);
   // contribution of bound constraints
   for (size_t i: this->lower_bounded_variables) {
      this->objective_gradient[i] -= this->barrier_parameter / (current_iterate.x[i] - this->variables_bounds[i].lb);
   }
   for (size_t i: this->upper_bounded_variables) {
      this->objective_gradient[i] -= this->barrier_parameter / (current_iterate.x[i] - this->variables_bounds[i].ub);
   }

   this->update_objective_multiplier(problem, current_iterate, objective_multiplier);

   /* bounds of the variables */
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->set_constraints_bounds(problem, this->constraints);

   /* set the initial point */
   clear(this->initial_point);
}

void InteriorPoint::update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) {
   // evaluate the Hessian
   this->hessian_evaluation->compute(problem, current_iterate.x, objective_multiplier, this->constraints_multipliers);

   // scale objective gradient
   if (objective_multiplier == 0.) {
      clear(this->objective_gradient);
   }
   else if (objective_multiplier < 1.) {
      scale(this->objective_gradient, objective_multiplier);
   }
}

/* reduced primal-dual approach */
Direction InteriorPoint::compute_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   this->iteration++;
   DEBUG << "\nCurrent iterate: " << current_iterate;

   // update the barrier parameter if the current iterate solves the subproblem
   this->update_barrier_parameter(problem, current_iterate);
   DEBUG << "mu is " << this->barrier_parameter << "\n";

   /************************/
   /* solve the KKT system */
   /************************/
   /* assemble and factorize the KKT matrix */
   COOMatrix kkt_matrix = this->assemble_kkt_matrix(problem, current_iterate);
   this->factorize(kkt_matrix, problem.type);

   /* inertia correction */
   this->modify_inertia(kkt_matrix, current_iterate.x.size(), problem.number_constraints, problem.type);
   DEBUG << "KKT matrix:\n" << kkt_matrix << "\n";

   /* right-hand side */
   this->generate_kkt_rhs(problem, current_iterate.x);
   std::vector<double> rhs_copy(rhs);

   /* compute the solution (Δx, -Δλ) */
   this->linear_solver->solve(this->rhs);
   this->number_subproblems_solved++;
   std::vector<double>& solution_IPM = this->rhs;

   /* generate IPM direction */
   Direction direction = this->generate_direction(problem, current_iterate, solution_IPM);
   direction.status = OPTIMAL;
   direction.norm = norm_inf(direction.x, 0, this->number_variables);
   /* evaluate the barrier objective */
   direction.objective = InteriorPoint::compute_barrier_directional_derivative(problem, current_iterate, direction.x);
   direction.predicted_reduction = [&](double step_length) {
      return this->compute_predicted_reduction(direction, step_length);
   };
   //std::cout << "IPM direction:\n" << direction << "\n";

   statistics.add_statistic("barrier param.", this->barrier_parameter);

   // check the residuals
//      std::vector<double> result = kkt_matrix.product(solution_IPM);
//      std::cout << "KKT matrix vector product: "; print_vector(std::cout, result);
//      std::cout << "RHS:                       "; print_vector(std::cout, rhs_copy);
//      std::vector<double> residual = add_vectors(result, rhs_copy, -1.);
//      double residual_norm = 0.;
//      for (size_t k = this->number_variables; k < residual.size(); k++) {
//         residual_norm += std::abs(residual[k]);
//      }
//      std::cout << "Residual norm: " << residual_norm << "\n";

   return direction;
//   catch (const UnstableInertiaCorrection& e) {
//      /* unstable factorization during optimality phase */
//      throw "InteriorPoint: inertia correction failed";
//   }
}

void InteriorPoint::update_barrier_parameter(const Problem& problem, Iterate& current_iterate) {
   /* scaled error terms */
   this->compute_residuals(problem, current_iterate, current_iterate.multipliers, 1.);
   double sd = this->compute_KKT_error_scaling(current_iterate);
   double KKTerror = current_iterate.residuals.KKT / sd;
   double central_complementarity_error = this->compute_central_complementarity_error(current_iterate, this->barrier_parameter, this->variables_bounds);
   DEBUG << "IPM error (KKT: " << KKTerror << ", cmpl: " << central_complementarity_error << ", feas: " << current_iterate.residuals.constraints <<
         ")\n";

   /* update of the barrier problem */
   double error = std::max(KKTerror, std::max(central_complementarity_error, current_iterate.residuals.constraints));
   if (error <= this->parameters.k_epsilon * this->barrier_parameter) {
      // TODO pass tolerance
      double tolerance = 1e-8;
      this->barrier_parameter = std::max(tolerance / 10.,
            std::min(this->parameters.k_mu * this->barrier_parameter, std::pow(this->barrier_parameter, this->parameters.theta_mu)));
      DEBUG << "IPM: mu updated to " << this->barrier_parameter << " and filter reset\n";
      // signal the redefinition of the problem to the globalization strategy
      this->subproblem_definition_changed = true;
   }
}

COOMatrix InteriorPoint::assemble_kkt_matrix(const Problem& problem, Iterate& current_iterate) {
   /* compute the Lagrangian Hessian */
   COOMatrix kkt_matrix = this->hessian_evaluation->hessian.to_COO();
   kkt_matrix.dimension = this->number_variables + problem.number_constraints;

   /* diagonal terms: bounds of primals and slacks */
   for (size_t i: this->lower_bounded_variables) {
      kkt_matrix.insert(current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - this->variables_bounds[i].lb), i, i);
   }
   for (size_t i: this->upper_bounded_variables) {
      kkt_matrix.insert(current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - this->variables_bounds[i].ub), i, i);
   }

   /* Jacobian of general constraints */
   for (size_t j = 0; j < problem.number_constraints; j++) {
      for (const auto[variable_index, derivative]: this->constraints_jacobian[j]) {
         kkt_matrix.insert(derivative, variable_index, this->number_variables + j);
      }
   }
   return kkt_matrix;
}

void InteriorPoint::generate_kkt_rhs(const Problem& problem, const std::vector<double>& current_primals) {
   /* generate the right-hand side */
   clear(this->rhs);

   /* barrier objective gradient */
   for (const auto[i, derivative]: this->objective_gradient) {
      this->rhs[i] = -derivative;
   }

   /* constraint gradients */
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (this->constraints_multipliers[j] != 0.) {
         for (const auto[variable_index, derivative]: this->constraints_jacobian[j]) {
            this->rhs[variable_index] += this->constraints_multipliers[j] * derivative;
         }
      }
   }

   /* constraint evaluations */
   int slack_index = 0;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (problem.constraint_status[j] == EQUAL_BOUNDS) { // add the bound
         this->rhs[this->number_variables + j] = problem.constraint_bounds[j].lb - constraints[j];
      }
      else { // add the slack
         this->rhs[this->number_variables + j] = current_primals[problem.number_variables + slack_index] - constraints[j];
         slack_index++;
      }
   }
   DEBUG << "RHS: ";
   print_vector(DEBUG, this->rhs);
}

double InteriorPoint::compute_KKT_error_scaling(Iterate& current_iterate) const {
   /* KKT error */
   double norm_1_constraint_multipliers = norm_1(current_iterate.multipliers.constraints);
   double norm_1_bound_multipliers = norm_1(current_iterate.multipliers.lower_bounds) + norm_1(current_iterate.multipliers.upper_bounds);
   double sd = std::max(this->parameters.smax, (norm_1_constraint_multipliers + norm_1_bound_multipliers) /
                                               (double) (current_iterate.x.size() + current_iterate.multipliers.constraints.size())) / this->parameters.smax;
   return sd;
}

Direction InteriorPoint::generate_direction(const Problem& problem, const Iterate& current_iterate, std::vector<double>& solution_IPM) {
   /* retrieve +Δλ (Nocedal p590) */
   for (size_t j = 0; j < problem.number_constraints; j++) {
      size_t multiplier_index = this->number_variables + j;
      solution_IPM[multiplier_index] = -solution_IPM[multiplier_index];
   }

   /* "fraction to boundary" rule for variables and bound multipliers */
   std::vector<double> trial_x(current_iterate.x.size());
   Multipliers trial_multipliers(current_iterate.x.size(), problem.number_constraints);
   double tau = std::max(this->parameters.tau_min, 1. - this->barrier_parameter);
   // scale primal variables and constraints multipliers
   double primal_length = this->compute_primal_length(current_iterate, solution_IPM, tau);
   for (size_t i = 0; i < this->number_variables; i++) {
      trial_x[i] = primal_length * solution_IPM[i];
   }
   for (size_t j = 0; j < problem.number_constraints; j++) {
      trial_multipliers.constraints[j] = current_iterate.multipliers.constraints[j] + primal_length * solution_IPM[number_variables + j];
   }

   /* compute bound multiplier displacements Δz */
   std::vector<double> lower_delta_z = this->compute_lower_bound_dual_displacements(current_iterate, solution_IPM);
   std::vector<double> upper_delta_z = this->compute_upper_bound_dual_displacements(current_iterate, solution_IPM);
   // scale dual variables
   double dual_length = InteriorPoint::compute_dual_length(current_iterate, tau, lower_delta_z, upper_delta_z);
   for (size_t i = 0; i < number_variables; i++) {
      trial_multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_length * lower_delta_z[i];
      trial_multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_length * upper_delta_z[i];
      // rescale the multipliers (IPOPT paper p6)
      //trial_multipliers.lower_bounds[i] = std::max(std::min(trial_multipliers.lower_bounds[i], this->parameters.kappa*this->mu_optimality/(trial_x[i] - subproblem_variables_bounds[i].lb)), this->mu_optimality/(this->parameters.kappa*(trial_x[i] - subproblem_variables_bounds[i].lb)));
   }

   DEBUG << "MA57 solution:\n";
   DEBUG << "Δx: ";
   print_vector(DEBUG, solution_IPM, '\n',0, problem.number_variables);
   DEBUG << "Δs: ";
   print_vector(DEBUG, solution_IPM, '\n', problem.number_variables, problem.inequality_constraints.size());
   DEBUG << "Δz_L: ";
   print_vector(DEBUG, lower_delta_z);
   DEBUG << "Δz_U: ";
   print_vector(DEBUG, upper_delta_z);
   DEBUG << "Δλ: ";
   print_vector(DEBUG, solution_IPM,'\n', number_variables, problem.number_constraints);
   DEBUG << "primal length = " << primal_length << "\n";
   DEBUG << "dual length = " << dual_length << "\n\n";

   Direction direction(trial_x, trial_multipliers);
   return direction;
}

double
InteriorPoint::compute_primal_length(const Iterate& current_iterate, const std::vector<double>& ipm_solution, double tau) {
   double primal_length = 1.;
   for (size_t i: this->lower_bounded_variables) {
      double trial_alpha_xi = -tau * (current_iterate.x[i] - this->variables_bounds[i].lb) / ipm_solution[i];
      if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   for (size_t i: this->upper_bounded_variables) {
      double trial_alpha_xi = -tau * (current_iterate.x[i] - this->variables_bounds[i].ub) / ipm_solution[i];
      if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   return primal_length;
}

double InteriorPoint::compute_dual_length(const Iterate& current_iterate, double tau, const std::vector<double>& lower_delta_z,
      const std::vector<double>& upper_delta_z) {
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

void InteriorPoint::factorize(COOMatrix& kkt_matrix, FunctionType problem_type) {
   // compute the symbolic factorization only when:
   // the problem has a non constant Hessian (ie is not an LP or a QP) or it is the first factorization
   // TODO: for QPs as well, but only when the sparsity pattern is constant
   if (force_symbolic_factorization || problem_type == LINEAR || this->number_factorizations_ == 0) {
      this->linear_solver->do_symbolic_factorization(kkt_matrix);
   }
   this->linear_solver->do_numerical_factorization(kkt_matrix);
   this->number_factorizations_++;
}

void InteriorPoint::modify_inertia(COOMatrix& kkt_matrix, size_t size_first_block, size_t size_second_block, FunctionType problem_type) {
   this->inertia_hessian = 0.;
   this->inertia_constraints = 0.;
   DEBUG << "Testing factorization with inertia term " << this->inertia_hessian << "\n";

   bool good_inertia = false;
   if (!this->linear_solver->matrix_is_singular() && this->linear_solver->number_negative_eigenvalues() == size_second_block) {
      DEBUG << "Factorization was a success\n";
      good_inertia = true;
   }
   else {
      // inertia term for constraints
      if (this->linear_solver->matrix_is_singular()) {
         DEBUG << "Matrix is singular\n";
         this->inertia_constraints = 1e-8 * std::pow(this->barrier_parameter, 0.25);
      }
      else {
         this->inertia_constraints = 0.;
      }
      // inertia term for Hessian
      if (this->inertia_hessian_last_ == 0.) {
         this->inertia_hessian = 1e-4;
      }
      else {
         this->inertia_hessian = std::max(1e-20, this->inertia_hessian_last_ / 3.);
      }
   }

   size_t current_matrix_size = kkt_matrix.matrix.size();
   if (!good_inertia) {
      for (size_t i = 0; i < size_first_block; i++) {
         kkt_matrix.insert(this->inertia_hessian, i, i);
      }
      for (size_t j = size_first_block; j < size_first_block + size_second_block; j++) {
         kkt_matrix.insert(-this->inertia_constraints, j, j);
      }
   }

   while (!good_inertia) {
      DEBUG << "Testing factorization with inertia term " << this->inertia_hessian << "\n";
      this->factorize(kkt_matrix, problem_type);

      if (!this->linear_solver->matrix_is_singular() && this->linear_solver->number_negative_eigenvalues() == size_second_block) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n";
         this->inertia_hessian_last_ = this->inertia_hessian;
      }
      else {
         if (this->inertia_hessian_last_ == 0.) {
            this->inertia_hessian *= 100.;
         }
         else {
            this->inertia_hessian *= 8.;
         }
         if (1e40 < this->inertia_hessian) {
            throw UnstableInertiaCorrection();
         }
         else {
            for (size_t i = 0; i < size_first_block; i++) {
               kkt_matrix.matrix[current_matrix_size + i] = this->inertia_hessian;
            }
            for (size_t j = size_first_block; j < size_first_block + size_second_block; j++) {
               kkt_matrix.matrix[current_matrix_size + j] = -this->inertia_constraints;
            }
         }
      }
   }
}

std::vector<double> InteriorPoint::compute_lower_bound_dual_displacements(const Iterate& current_iterate, const std::vector<double>& solution) {
   std::vector<double> delta_z(current_iterate.multipliers.lower_bounds.size());
   for (size_t i: this->lower_bounded_variables) {
      delta_z[i] = this->barrier_parameter / (current_iterate.x[i] - this->variables_bounds[i].lb) - current_iterate.multipliers.lower_bounds[i] -
                   current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - this->variables_bounds[i].lb) * solution[i];
   }
   return delta_z;
}

std::vector<double> InteriorPoint::compute_upper_bound_dual_displacements(const Iterate& current_iterate, const std::vector<double>& solution) {
   std::vector<double> delta_z(current_iterate.multipliers.upper_bounds.size());
   for (size_t i: this->upper_bounded_variables) {
      delta_z[i] = this->barrier_parameter / (current_iterate.x[i] - this->variables_bounds[i].ub) - current_iterate.multipliers.upper_bounds[i] -
                   current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - this->variables_bounds[i].ub) * solution[i];
   }
   return delta_z;
}

void InteriorPoint::compute_progress_measures(const Problem& problem, Iterate& iterate) {
   /* evaluate constraints with slacks */
   double feasibility = InteriorPoint::constraint_violation(problem, iterate);
   /* compute barrier objective */
   double objective = this->evaluate_barrier_function(problem, iterate);
   iterate.progress = {feasibility, objective};
}

double InteriorPoint::constraint_violation(const Problem& problem, Iterate& iterate) {
   iterate.compute_constraints(problem);
   // compute norm by using a lambda (avoids constructing a vector)
   size_t slack_index = problem.number_variables;
   auto residual_function = [&](size_t j) {
      if (problem.constraint_status[j] == EQUAL_BOUNDS) {
         return iterate.constraints[j] - problem.constraint_bounds[j].lb;
      }
      else {
         double constraint_value = iterate.constraints[j] - iterate.x[slack_index];
         slack_index++;
         return constraint_value;
      }
   };
   return norm(residual_function, problem.number_constraints, L1_NORM);
}

double InteriorPoint::evaluate_barrier_function(const Problem& problem, Iterate& iterate) {
   double objective = 0.;
   /* bound constraints */
   for (size_t i: this->lower_bounded_variables) {
      objective -= std::log(iterate.x[i] - this->variables_bounds[i].lb);
   }
   for (size_t i: this->upper_bounded_variables) {
      objective -= std::log(this->variables_bounds[i].ub - iterate.x[i]);
   }
   objective *= this->barrier_parameter;
   /* original objective */
   iterate.compute_objective(problem);
   objective += iterate.objective;
   return objective;
}

double InteriorPoint::compute_barrier_directional_derivative(const Problem& /*problem*/, Iterate& current_iterate, std::vector<double>& solution) {
   double subproblem_objective = dot(solution, current_iterate.objective_gradient);
   return subproblem_objective;
}

double InteriorPoint::compute_predicted_reduction(const Direction& direction, double step_length) const {
   // the predicted reduction is linear
   return -step_length * direction.objective;
}

double InteriorPoint::compute_central_complementarity_error(Iterate& iterate, double mu, std::vector<Range>& variables_bounds) const {
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
   double sc = std::max(this->parameters.smax,
         (norm_1(iterate.multipliers.lower_bounds) + norm_1(iterate.multipliers.upper_bounds)) / (double) iterate.x.size()) / this->parameters.smax;
   return norm(residuals, L1_NORM) / sc;
}

int InteriorPoint::get_hessian_evaluation_count() const {
   return this->hessian_evaluation->evaluation_count;
}