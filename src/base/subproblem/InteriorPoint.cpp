#include <cmath>
#include <cassert>
#include <iomanip>
#include "InteriorPoint.hpp"
#include "LinearSolverFactory.hpp"

InteriorPoint::InteriorPoint(const Problem& problem, size_t number_variables, size_t number_constraints, const std::string& linear_solver_name, const
std::string& hessian_evaluation_method, bool use_trust_region) :
      // add the slacks to the variables
      Subproblem(number_variables + problem.inequality_constraints.size(), number_constraints),
      /* if no trust region is used, the problem should be convexified. However, the inertia of the augmented matrix will be corrected later */
      hessian_evaluation(HessianEvaluationFactory::create(hessian_evaluation_method, problem.number_variables, problem
            .hessian_maximum_number_nonzeros, false)),
      kkt_matrix(this->number_variables + number_constraints, problem.hessian_maximum_number_nonzeros, 1),
      linear_solver(LinearSolverFactory::create(linear_solver_name)),
      parameters({0.99, 1e10, 100., 0.2, 1.5, 10., 1e10}),
      rhs(this->number_variables + number_constraints),
      lower_delta_z(this->number_variables), upper_delta_z(this->number_variables) {
   // register the original variables bounds
   for (size_t i = 0; i < problem.number_variables; i++) {
      this->variables_bounds[i] = problem.variables_bounds[i];
   }

   // identify the bounded variables
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
      size_t slack_index = number_variables + i;
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

void InteriorPoint::evaluate_constraints(const Problem& problem, Iterate& iterate) const {
   iterate.compute_constraints(problem);
   // transform the constraints into "= 0" equalities
   for (const auto& element: problem.equality_constraints) {
      size_t j = element.first;
      iterate.constraints[j] -= problem.constraint_bounds[j].lb;
   }
   for (const auto[j, i]: problem.inequality_constraints) {
      iterate.constraints[j] -= iterate.x[problem.number_variables + i];
   }
}

Iterate InteriorPoint::generate_initial_iterate(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
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
      multipliers.lower_bounds[i] = this->default_multiplier_;
   }
   for (size_t i: this->upper_bounded_variables) {
      multipliers.upper_bounds[i] = -this->default_multiplier_;
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
   this->evaluate_constraints(problem, first_iterate);
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

void InteriorPoint::generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);
   /* compute first- and second-order information */
   // constraint Jacobian
   for (auto& row: this->constraints_jacobian) {
      row.clear();
   }
   problem.constraints_jacobian(current_iterate.x, this->constraints_jacobian);
   // add the slack variables
   for (const auto[j, i]: problem.inequality_constraints) {
      this->constraints_jacobian[j][problem.number_variables + i] = -1.;
   }

   // objective gradient
   this->objective_gradient.clear();
   problem.objective_gradient(current_iterate.x, this->objective_gradient);
   for (size_t i: this->lower_bounded_variables) {
      this->objective_gradient[i] -= this->barrier_parameter / (current_iterate.x[i] - this->variables_bounds[i].lb);
   }
   for (size_t i: this->upper_bounded_variables) {
      this->objective_gradient[i] -= this->barrier_parameter / (current_iterate.x[i] - this->variables_bounds[i].ub);
   }

   // Hessian (scaled by the objective multiplier)
   this->update_objective_multiplier(problem, current_iterate, objective_multiplier);

   /* bounds of the variables */
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);
}

void InteriorPoint::update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) {
   // evaluate the Hessian
   this->hessian_evaluation->compute(problem, current_iterate.x, objective_multiplier, this->constraints_multipliers);

   // scale objective gradient
   if (objective_multiplier == 0.) {
      clear(this->objective_gradient);
   }
   else if (objective_multiplier < 1.) {
      this->objective_gradient = current_iterate.objective_gradient;
      scale(this->objective_gradient, objective_multiplier);
   }
}

void InteriorPoint::set_initial_point(const std::vector<double>& /*point*/) {
   // do nothing
}

/* reduced primal-dual approach */
Direction InteriorPoint::compute_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   this->iteration++;
   // update the barrier parameter if the current iterate solves the subproblem
   this->update_barrier_parameter(current_iterate);
   DEBUG << "mu is " << this->barrier_parameter << "\n";

   /************************/
   /* solve the KKT system */
   /************************/
   /* assemble and factorize the KKT matrix */
   this->kkt_matrix = this->assemble_kkt_matrix(problem, current_iterate);
   this->factorize(this->kkt_matrix, problem.type);
   /* inertia correction */
   this->modify_inertia(this->kkt_matrix, this->number_variables, problem.number_constraints, problem.type);
   DEBUG << "KKT matrix:\n" << this->kkt_matrix << "\n";
   auto [number_pos, number_neg, number_zero] = this->linear_solver->get_inertia();
   assert(number_pos == this->number_variables && number_neg == problem.number_constraints && number_zero == 0);

   /* right-hand side */
   this->generate_kkt_rhs(current_iterate);

   /* compute the solution (Δx, -Δλ) */
   std::vector<double> solution_IPM = this->linear_solver->solve(this->kkt_matrix, this->rhs);
   this->number_subproblems_solved++;

   /* generate IPM direction */
   Direction direction = this->generate_direction(problem, current_iterate, solution_IPM);
   direction.status = OPTIMAL;
   direction.norm = norm_inf(direction.x, 0, this->number_variables);
   /* evaluate the barrier objective */
   direction.objective = this->compute_barrier_directional_derivative(direction.x);
   //assert(direction.objective < 0. && "the IPM directional derivative is positive");

   statistics.add_statistic("barrier param.", this->barrier_parameter);
   return direction;
   //   catch (const UnstableInertiaCorrection& e) {
   //      /* unstable factorization during optimality phase */
   //      throw "InteriorPoint: inertia correction failed";
   //   }
}

void InteriorPoint::update_barrier_parameter(const Iterate& current_iterate) {
   /* scaled error terms */
   double sd = this->compute_KKT_error_scaling(current_iterate);
   double KKTerror = current_iterate.errors.KKT / sd;
   double central_complementarity_error = this->compute_central_complementarity_error(current_iterate);
   DEBUG << "IPM error (KKT: " << KKTerror << ", cmpl: " << central_complementarity_error << ", feas: " << current_iterate.errors.constraints <<
         ")\n";

   /* update of the barrier problem */
   double error = std::max(KKTerror, std::max(central_complementarity_error, current_iterate.errors.constraints));
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
      for (const auto[i, derivative]: this->constraints_jacobian[j]) {
         kkt_matrix.insert(derivative, i, this->number_variables + j);
      }
   }
   return kkt_matrix;
}

void InteriorPoint::generate_kkt_rhs(const Iterate& current_iterate) {
   /* generate the right-hand side */
   clear(this->rhs);

   /* barrier objective gradient */
   for (const auto[i, derivative]: this->objective_gradient) {
      this->rhs[i] = -derivative;
   }

   /* constraint: evaluations and gradients */
   for (size_t j = 0; j < current_iterate.constraints.size(); j++) {
      // Lagrangian
      if (this->constraints_multipliers[j] != 0.) {
         for (const auto[i, derivative]: this->constraints_jacobian[j]) {
            this->rhs[i] += this->constraints_multipliers[j] * derivative;
         }
      }
      // constraints
      this->rhs[this->number_variables + j] = -current_iterate.constraints[j];
   }
   DEBUG << "RHS: "; print_vector(DEBUG, this->rhs); DEBUG << "\n";
}

Direction InteriorPoint::compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) {
   DEBUG << "Entered SOC computation\n";
   DEBUG << "KKT matrix:\n" << this->kkt_matrix << "\n";
   for (size_t j = 0; j < trial_iterate.constraints.size(); j++) {
      this->rhs[this->number_variables + j] -= trial_iterate.constraints[j];
   }
   DEBUG << "SOC RHS: "; print_vector(DEBUG, this->rhs); DEBUG << "\n";

   /* compute the solution (Δx, -Δλ) */
   std::vector<double> solution_IPM = this->linear_solver->solve(kkt_matrix, this->rhs);
   this->number_subproblems_solved++;

   /* generate IPM direction */
   Direction direction_soc = this->generate_direction(problem, trial_iterate, solution_IPM);
   direction_soc.status = OPTIMAL;
   direction_soc.norm = norm_inf(direction_soc.x, 0, this->number_variables);
   /* evaluate the barrier objective */
   direction_soc.objective = this->compute_barrier_directional_derivative(direction_soc.x);
   DEBUG << "SOC direction:\n" << direction_soc << "\n";
   return direction_soc;
}

double InteriorPoint::compute_KKT_error_scaling(const Iterate& current_iterate) const {
   /* KKT error */
   const double norm_1_constraint_multipliers = norm_1(current_iterate.multipliers.constraints);
   const double norm_1_bound_multipliers = norm_1(current_iterate.multipliers.lower_bounds) + norm_1(current_iterate.multipliers.upper_bounds);
   const double sd = std::max(this->parameters.smax, (norm_1_constraint_multipliers + norm_1_bound_multipliers) /
       (double) (current_iterate.x.size() + current_iterate.multipliers.constraints.size())) / this->parameters.smax;
   return sd;
}

Direction InteriorPoint::generate_direction(const Problem& problem, const Iterate& current_iterate, std::vector<double>& solution_IPM) {
   /* retrieve +Δλ (Nocedal p590) */
   for (size_t j = this->number_variables; j < solution_IPM.size(); j++) {
      solution_IPM[j] = -solution_IPM[j];
   }

   Direction direction(this->number_variables, problem.number_constraints);

   // "fraction to boundary" rule for primal variables and constraints multipliers
   double tau = this->parameters.tau_min; //std::max(this->parameters.tau_min, 1. - this->barrier_parameter);
   double primal_step_length = this->primal_fraction_to_boundary(current_iterate, solution_IPM, tau);
   for (size_t i = 0; i < this->number_variables; i++) {
      direction.x[i] = primal_step_length * solution_IPM[i];
   }
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] = primal_step_length * solution_IPM[number_variables + j];
   }

   /* compute bound multiplier displacements Δz */
   this->compute_lower_bound_dual_displacements(current_iterate, solution_IPM);
   this->compute_upper_bound_dual_displacements(current_iterate, solution_IPM);

   // "fraction to boundary" rule for bound multipliers
   double dual_step_length = this->dual_fraction_to_boundary(current_iterate, tau);
   for (size_t i = 0; i < number_variables; i++) {
      direction.multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_step_length * this->lower_delta_z[i];
      direction.multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_step_length * this->upper_delta_z[i];
   }

   DEBUG << "IPM solution:\n";
   DEBUG << "Δx: "; print_vector(DEBUG, solution_IPM, '\n', 0, problem.number_variables);
   DEBUG << "Δs: "; print_vector(DEBUG, solution_IPM, '\n', problem.number_variables, problem.inequality_constraints.size());
   DEBUG << "Δz_L: "; print_vector(DEBUG, this->lower_delta_z);
   DEBUG << "Δz_U: "; print_vector(DEBUG, this->upper_delta_z);
   DEBUG << "Δλ: "; print_vector(DEBUG, solution_IPM, '\n', number_variables, problem.number_constraints);
   DEBUG << "primal length = " << primal_step_length << "\n";
   DEBUG << "dual length = " << dual_step_length << "\n\n";
   return direction;
}

double
InteriorPoint::primal_fraction_to_boundary(const Iterate& current_iterate, const std::vector<double>& ipm_solution, double tau) {
   double primal_length = 1.;
   for (size_t i: this->lower_bounded_variables) {
      if (ipm_solution[i] < 0.) {
         double trial_alpha_xi = -tau * (current_iterate.x[i] - this->variables_bounds[i].lb) / ipm_solution[i];
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   for (size_t i: this->upper_bounded_variables) {
      if (0. < ipm_solution[i]) {
         double trial_alpha_xi = -tau * (current_iterate.x[i] - this->variables_bounds[i].ub) / ipm_solution[i];
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   return primal_length;
}

double InteriorPoint::dual_fraction_to_boundary(const Iterate& current_iterate, double tau) {
   double dual_length = 1.;
   for (size_t i = 0; i < current_iterate.multipliers.lower_bounds.size(); i++) {
      if (this->lower_delta_z[i] < 0.) {
         double trial_alpha_zj = -tau * current_iterate.multipliers.lower_bounds[i] / this->lower_delta_z[i];
         dual_length = std::min(dual_length, trial_alpha_zj);
      }
      if (0. < this->upper_delta_z[i]) {
         double trial_alpha_zj = -tau * current_iterate.multipliers.upper_bounds[i] / this->upper_delta_z[i];
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

void InteriorPoint::compute_lower_bound_dual_displacements(const Iterate& current_iterate, const std::vector<double>& solution) {
   clear(this->lower_delta_z);
   for (size_t i: this->lower_bounded_variables) {
      this->lower_delta_z[i] = (this->barrier_parameter - solution[i]*current_iterate.multipliers.lower_bounds[i]) / (current_iterate.x[i] -
            this->variables_bounds[i].lb) - current_iterate.multipliers.lower_bounds[i];
   }
}

void InteriorPoint::compute_upper_bound_dual_displacements(const Iterate& current_iterate, const std::vector<double>& solution) {
   clear(this->upper_delta_z);
   for (size_t i: this->upper_bounded_variables) {
      this->upper_delta_z[i] = (this->barrier_parameter - solution[i]*current_iterate.multipliers.upper_bounds[i]) / (current_iterate.x[i] -
            this->variables_bounds[i].ub) - current_iterate.multipliers.upper_bounds[i];
   }
}

double InteriorPoint::compute_constraint_violation(const Problem& /*problem*/, const Iterate& iterate) const {
   return norm_1(iterate.constraints);
}

void InteriorPoint::compute_progress_measures(const Problem& problem, Iterate& iterate) {
   const double constraint_violation = this->compute_constraint_violation(problem, iterate);
   /* compute barrier objective */
   const double objective = this->evaluate_barrier_function(problem, iterate);
   iterate.progress = {constraint_violation, objective};
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

double InteriorPoint::compute_barrier_directional_derivative(const std::vector<double>& solution) {
   return dot(solution, this->objective_gradient);
}

double InteriorPoint::compute_predicted_reduction(const Direction& direction, double step_length) const {
   // the predicted reduction is linear
   return -step_length * direction.objective;
}

double InteriorPoint::compute_central_complementarity_error(const Iterate& iterate) const {
   /* variable bound constraints */
   const auto residual_function = [&](size_t i) {
      double result = 0.;
      if (-INFINITY < this->variables_bounds[i].lb) {
         result += iterate.multipliers.lower_bounds[i] * (iterate.x[i] - this->variables_bounds[i].lb) - this->barrier_parameter;
      }
      if (this->variables_bounds[i].ub < INFINITY) {
         result += iterate.multipliers.upper_bounds[i] * (iterate.x[i] - this->variables_bounds[i].ub) - this->barrier_parameter;
      }
      return result;
   };

   /* scaling */
   double sc = std::max(this->parameters.smax,
         (norm_1(iterate.multipliers.lower_bounds) + norm_1(iterate.multipliers.upper_bounds)) / (double) iterate.x.size()) / this->parameters.smax;
   return norm_1(residual_function, iterate.x.size()) / sc;
}

int InteriorPoint::get_hessian_evaluation_count() const {
   return this->hessian_evaluation->evaluation_count;
}
