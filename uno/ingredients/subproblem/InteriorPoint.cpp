#include "InteriorPoint.hpp"
#include "solvers/linear/LinearSolverFactory.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"
#include "optimization/Preprocessing.hpp"

InteriorPoint::InteriorPoint(const Problem& problem, size_t max_number_variables, size_t number_constraints,
      const std::string& hessian_model, const std::string& linear_solver_name, const std::string& sparse_format, double initial_barrier_parameter,
      double default_multiplier, double tolerance, bool use_trust_region, const Options& options) :
      // add the slacks to the variables
      Subproblem(problem.number_variables + problem.inequality_constraints.size(), // number_variables
            max_number_variables + problem.inequality_constraints.size(), // max_number_variables
            number_constraints, SOC_UPON_REJECTION, true),
      augmented_system(sparse_format, this->max_number_variables + number_constraints,
            problem.hessian_maximum_number_nonzeros
            + this->max_number_variables + number_constraints /* regularization */
            + 2 * this->max_number_variables /* diagonal barrier terms */
            + this->max_number_variables * number_constraints /* Jacobian */,
            stod(options.at("LS_regularization_failure_threshold"))),
      barrier_parameter(initial_barrier_parameter), tolerance(tolerance),
      // if no trust region is used, the problem should be convexified. However, the inertia of the augmented matrix will be corrected later
      hessian_model(HessianModelFactory::create(hessian_model, this->max_number_variables, problem.hessian_maximum_number_nonzeros, sparse_format,
            false)),
      linear_solver(LinearSolverFactory::create(linear_solver_name, this->max_number_variables + number_constraints,
            problem.hessian_maximum_number_nonzeros
            + this->max_number_variables + number_constraints /* regularization */
            + 2 * this->max_number_variables /* diagonal barrier terms */
            + this->max_number_variables * number_constraints /* Jacobian */)),
      parameters({stod(options.at("tau_min")),
            stod(options.at("k_sigma")),
            stod(options.at("smax")),
            stod(options.at("k_mu")),
            stod(options.at("theta_mu")),
            stod(options.at("k_epsilon")),
            stod(options.at("barrier_update_fraction")),
            stod(options.at("regularization_barrier_exponent"))}),
      default_multiplier(default_multiplier),
      primal_iterate(this->max_number_variables),
      lower_bound_multipliers(this->max_number_variables),
      upper_bound_multipliers(this->max_number_variables),
      barrier_constraints(number_constraints),
      lower_delta_z(this->max_number_variables), upper_delta_z(this->max_number_variables) {
   // register the original variables bounds
   copy_from(this->variables_bounds, problem.variables_bounds, problem.number_variables);

   // constraints are transformed into "c(x) = 0"
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->constraints_bounds[j] = {0., 0.};
   }

   // identify the bounded variables
   for (size_t i = 0; i < problem.number_variables; i++) {
      if (use_trust_region || (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
         this->lower_bounded_variables.push_back(i);
      }
      if (use_trust_region || (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
         this->upper_bounded_variables.push_back(i);
      }
   }
   // identify the inequality constraint slacks
   DEBUG << problem.inequality_constraints.size() << " slacks\n";
   for (const auto[j, i]: problem.inequality_constraints) {
      const size_t slack_index = problem.number_variables + i;
      if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
         this->lower_bounded_variables.push_back(slack_index);
      }
      if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
         this->upper_bounded_variables.push_back(slack_index);
      }
      // store the bounds of the slacks
      this->variables_bounds[slack_index] = problem.constraint_bounds[j];
   }
}

void InteriorPoint::set_initial_point(const std::vector<double>& /*initial_point*/) {
   // do nothing
}

void InteriorPoint::set_constraints(const Problem& problem, Iterate& iterate) {
   iterate.evaluate_constraints(problem);
   // transform the constraints into "= 0" equalities
   for (const auto& element: problem.equality_constraints) {
      const size_t j = element.first;
      this->barrier_constraints[j] = iterate.constraints[j] - problem.constraint_bounds[j].lb;
   }
   for (const auto[j, i]: problem.inequality_constraints) {
      this->barrier_constraints[j] = iterate.constraints[j] - iterate.x[problem.number_variables + i];
   }
}

inline void InteriorPoint::initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) {
   statistics.add_column("barrier param.", Statistics::double_width, 8);

   // resize to the new size (primals + slacks)
   first_iterate.change_number_variables(this->number_variables);

   // make the initial point strictly feasible wrt the bounds
   for (size_t i = 0; i < problem.number_variables; i++) {
      first_iterate.x[i] = Subproblem::push_variable_to_interior(first_iterate.x[i], problem.variables_bounds[i]);
   }

   // initialize the slacks and add contribution to the constraint Jacobian
   first_iterate.evaluate_constraints(problem);
   first_iterate.evaluate_constraints_jacobian(problem);
   for (const auto[j, i]: problem.inequality_constraints) {
      const double slack_value = Subproblem::push_variable_to_interior(first_iterate.constraints[j], problem.constraint_bounds[j]);
      first_iterate.x[problem.number_variables + i] = slack_value;
      first_iterate.constraints_jacobian[j].insert(problem.number_variables + i, -1.);
   }
   this->set_current_iterate(first_iterate);

   // compute least-square multipliers
   if (problem.is_constrained()) {
      Preprocessing::compute_least_square_multipliers(problem, *this->augmented_system.matrix, this->augmented_system.rhs, *this->linear_solver,
            first_iterate, first_iterate.multipliers.constraints);
   }

   // set the bound multipliers
   for (size_t i: this->lower_bounded_variables) {
      first_iterate.multipliers.lower_bounds[i] = this->default_multiplier;
   }
   for (size_t i: this->upper_bounded_variables) {
      first_iterate.multipliers.upper_bounds[i] = -this->default_multiplier;
   }

   // compute the optimality and feasibility measures of the initial point
   this->set_constraints(problem, first_iterate);
   this->compute_progress_measures(problem, first_iterate);
}

void InteriorPoint::create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier,
      double trust_region_radius) {
   // update the barrier parameter if the current iterate solves the subproblem
   this->update_barrier_parameter(current_iterate);

   // save the current iterate locally
   this->set_current_iterate(current_iterate);

   // constraints
   this->set_constraints(problem, current_iterate);
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);

   // constraint Jacobian
   problem.evaluate_constraint_jacobian(current_iterate.x, this->constraints_jacobian);
   // add the slack variables
   for (const auto[j, i]: problem.inequality_constraints) {
      this->constraints_jacobian[j].insert(problem.number_variables + i, -1.);
   }

   // build a model of the objective scaled by the objective multiplier
   this->build_objective_model(problem, current_iterate, objective_multiplier);

   // variables and bounds
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);
}

void InteriorPoint::build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // evaluate the Hessian
   this->hessian_model->evaluate(problem, current_iterate.x, objective_multiplier, this->constraints_multipliers);

   // objective gradient
   this->set_scaled_objective_gradient(problem, current_iterate, objective_multiplier);
   for (size_t i: this->lower_bounded_variables) {
      this->objective_gradient.insert(i, -this->barrier_parameter / (this->primal_iterate[i] - this->variables_bounds[i].lb));
   }
   for (size_t i: this->upper_bounded_variables) {
      this->objective_gradient.insert(i, -this->barrier_parameter / (this->primal_iterate[i] - this->variables_bounds[i].ub));
   }
}

Direction InteriorPoint::solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   this->iteration++;
   // set up the augmented system (with the right inertia)
   this->hessian_model->finalize(this->number_variables);
   this->assemble_augmented_system(problem);

   // compute the solution (Δx, -Δλ)
   this->augmented_system.solve(*this->linear_solver, this->number_variables + this->number_constraints);
   this->number_subproblems_solved++;

   // generate IPM direction
   this->generate_direction(problem, current_iterate);
   statistics.add_statistic("barrier param.", this->barrier_parameter);
   return this->direction;
}

//   catch (const UnstableInertiaCorrection& e) {
//      // unstable factorization during optimality phase
//      throw "InteriorPoint: inertia correction failed";
//   }

void InteriorPoint::assemble_augmented_system(const Problem& problem) {
   // assemble, factorize and regularize the KKT matrix
   this->assemble_augmented_matrix();
   this->augmented_system.factorize_matrix(problem, *this->linear_solver, this->number_variables + this->number_constraints);
   this->augmented_system.regularize_matrix(problem, *this->linear_solver, this->number_variables, this->number_constraints,
         std::pow(this->barrier_parameter, this->parameters.regularization_barrier_exponent));
   auto[number_pos, number_neg, number_zero] = this->linear_solver->get_inertia();
   assert(number_pos == this->number_variables && number_neg == this->number_constraints && number_zero == 0);

   // right-hand side
   this->generate_augmented_rhs();
}

Direction InteriorPoint::compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) {
   // save the current iterate locally
   this->set_current_iterate(trial_iterate);

   // modify the RHS by adding the values of the constraints
   for (const auto& element: problem.equality_constraints) {
      size_t j = element.first;
      this->augmented_system.rhs[this->number_variables + j] -= trial_iterate.constraints[j] - problem.constraint_bounds[j].lb;
   }
   for (const auto[j, i]: problem.inequality_constraints) {
      this->augmented_system.rhs[this->number_variables + j] -= trial_iterate.constraints[j] - trial_iterate.x[problem.number_variables + i];
   }

   // compute the solution (Δx, -Δλ)
   this->augmented_system.solve(*this->linear_solver, this->number_variables + this->number_constraints);
   this->number_subproblems_solved++;

   // generate IPM direction
   this->generate_direction(problem, trial_iterate);
   this->print_soc_iteration(this->direction);
   return this->direction;
}

void InteriorPoint::add_variable(size_t i, double current_value, const Range& bounds, double objective_term, size_t j, double jacobian_term) {
   // add the variable to the objective and the constraint Jacobian
   Subproblem::add_variable(i, current_value, bounds, objective_term, j, jacobian_term);
   // if necessary, register the variable as bounded
   if (is_finite_lower_bound(bounds.lb)) {
      this->lower_bounded_variables.push_back(i);
      this->lower_bound_multipliers[i] = this->default_multiplier;
   }
   if (is_finite_upper_bound(bounds.ub)) {
      this->upper_bounded_variables.push_back(i);
      this->upper_bound_multipliers[i] = -this->default_multiplier;
   }
   // save the current value
   this->primal_iterate[i] = Subproblem::push_variable_to_interior(current_value, bounds);
}

void InteriorPoint::remove_variable(size_t i, size_t j) {
   // remove the variable to the objective and the constraint Jacobian
   Subproblem::remove_variable(i, j);
   this->lower_bounded_variables.erase(std::remove(this->lower_bounded_variables.begin(), this->lower_bounded_variables.end(), i),
         this->lower_bounded_variables.end());
   this->upper_bounded_variables.erase(std::remove(this->upper_bounded_variables.begin(), this->upper_bounded_variables.end(), i),
         this->upper_bounded_variables.end());
}

void InteriorPoint::print_soc_iteration(const Direction& direction_soc) const {
   DEBUG << "Entered SOC computation\n";
   DEBUG << "KKT matrix:\n" << *this->augmented_system.matrix << "\n";
   DEBUG << "SOC RHS: ";
   print_vector(DEBUG, this->augmented_system.rhs);
   DEBUG << "\nSOC direction:\n" << direction_soc << "\n";
}

PredictedReductionModel InteriorPoint::generate_predicted_reduction_model(const Problem& /*problem*/,
      const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() {
      return [=](double step_length) {
         return -step_length * direction.objective;
      };
   });
}

void InteriorPoint::compute_progress_measures(const Problem& problem, Iterate& iterate) {
   iterate.evaluate_constraints(problem);
   auto residual_function = [&](size_t j) {
      if (problem.constraint_status[j] == EQUAL_BOUNDS) {
         return iterate.constraints[j] - problem.constraint_bounds[j].lb;
      }
      else {
         const size_t i = problem.inequality_constraints.at(j);
         return iterate.constraints[j] - iterate.x[problem.number_variables + i];
      }
   };
   const double constraint_violation = norm_1(residual_function, problem.number_constraints);

   // compute barrier objective
   const double barrier_objective = this->evaluate_barrier_function(problem, iterate);
   iterate.progress = {constraint_violation, barrier_objective};
}

void InteriorPoint::update_barrier_parameter(const Iterate& current_iterate) {
   const double tolerance_fraction = this->tolerance / this->parameters.barrier_update_fraction;
   // scaled error terms
   const double sd = this->compute_KKT_error_scaling(current_iterate);
   const double KKTerror = current_iterate.errors.KKT / sd;
   const double central_complementarity_error = this->compute_central_complementarity_error(current_iterate);
   const double error = std::max({KKTerror, current_iterate.errors.constraints, central_complementarity_error});

   // update of the barrier problem
   while (error <= this->parameters.k_epsilon * this->barrier_parameter && tolerance_fraction < this->barrier_parameter) {
      this->barrier_parameter = std::max(tolerance_fraction, std::min(this->parameters.k_mu * this->barrier_parameter,
            std::pow(this->barrier_parameter, this->parameters.theta_mu)));
      DEBUG << "IPM: mu updated to " << this->barrier_parameter << " and filter reset\n";
      // signal the redefinition of the problem to the globalization strategy
      this->subproblem_definition_changed = true;
   }
   DEBUG << "mu is " << this->barrier_parameter << "\n";
}

void InteriorPoint::set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   // here, we work with the original bounds
   // very important: apply the trust region only on the original variables (not the slacks)
   for (size_t i = 0; i < problem.number_variables; i++) {
      double lb = std::max(current_iterate.x[i] - trust_region_radius, problem.variables_bounds[i].lb);
      double ub = std::min(current_iterate.x[i] + trust_region_radius, problem.variables_bounds[i].ub);
      this->variables_bounds[i] = {lb, ub};
   }
}

double InteriorPoint::compute_barrier_directional_derivative(const std::vector<double>& solution) {
   return dot(solution, this->objective_gradient);
}

double InteriorPoint::evaluate_barrier_function(const Problem& problem, Iterate& iterate) {
   double objective = 0.;
   // bound constraints
   for (size_t i: this->lower_bounded_variables) {
      objective -= std::log(iterate.x[i] - this->variables_bounds[i].lb);
   }
   for (size_t i: this->upper_bounded_variables) {
      objective -= std::log(this->variables_bounds[i].ub - iterate.x[i]);
   }
   objective *= this->barrier_parameter;
   // original objective
   iterate.evaluate_objective(problem);
   objective += iterate.objective;
   return objective;
}

double InteriorPoint::primal_fraction_to_boundary(const std::vector<double>& ipm_solution, double tau) {
   double primal_length = 1.;
   for (size_t i: this->lower_bounded_variables) {
      if (ipm_solution[i] < 0.) {
         double trial_alpha_xi = -tau * (this->primal_iterate[i] - this->variables_bounds[i].lb) / ipm_solution[i];
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   for (size_t i: this->upper_bounded_variables) {
      if (0. < ipm_solution[i]) {
         double trial_alpha_xi = -tau * (this->primal_iterate[i] - this->variables_bounds[i].ub) / ipm_solution[i];
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   return primal_length;
}

double InteriorPoint::dual_fraction_to_boundary(double tau) {
   double dual_length = 1.;
   for (size_t i = 0; i < this->number_variables; i++) {
      if (this->lower_delta_z[i] < 0.) {
         double trial_alpha_zj = -tau * this->lower_bound_multipliers[i] / this->lower_delta_z[i];
         dual_length = std::min(dual_length, trial_alpha_zj);
      }
      if (0. < this->upper_delta_z[i]) {
         double trial_alpha_zj = -tau * this->upper_bound_multipliers[i] / this->upper_delta_z[i];
         dual_length = std::min(dual_length, trial_alpha_zj);
      }
   }
   return dual_length;
}

void InteriorPoint::assemble_augmented_matrix() {
   this->augmented_system.matrix->reset();
   this->augmented_system.matrix->dimension = this->number_variables + this->number_constraints;
   // copy the Lagrangian Hessian in the top left block
   size_t current_column = 0;
   this->hessian_model->hessian->for_each([&](size_t i, size_t j, double entry) {
      for (size_t column = current_column; column < j; column++) {
         this->augmented_system.matrix->finalize(column);
         current_column++;
      }
      this->augmented_system.matrix->insert(entry, i, j);
   });

   // diagonal terms: bounds of primals and slacks
   for (size_t i: this->lower_bounded_variables) {
      this->augmented_system.matrix->insert(this->lower_bound_multipliers[i] / (this->primal_iterate[i] - this->variables_bounds[i].lb), i, i);
   }
   for (size_t i: this->upper_bounded_variables) {
      this->augmented_system.matrix->insert(this->upper_bound_multipliers[i] / (this->primal_iterate[i] - this->variables_bounds[i].ub), i, i);
   }

   // Jacobian of general constraints
   for (size_t j = 0; j < this->number_constraints; j++) {
      this->constraints_jacobian[j].for_each([&](size_t i, double derivative) {
         this->augmented_system.matrix->insert(derivative, i, this->number_variables + j);
      });
      this->augmented_system.matrix->finalize(j);
   }
}

void InteriorPoint::generate_augmented_rhs() {
   // generate the right-hand side
   clear(this->augmented_system.rhs);

   // barrier objective gradient
   this->objective_gradient.for_each([&](size_t i, double derivative) {
      this->augmented_system.rhs[i] = -derivative;
   });

   // constraint: evaluations and gradients
   for (size_t j = 0; j < this->number_constraints; j++) {
      // Lagrangian
      if (this->constraints_multipliers[j] != 0.) {
         this->constraints_jacobian[j].for_each([&](size_t i, double derivative) {
            this->augmented_system.rhs[i] += this->constraints_multipliers[j] * derivative;
         });
      }
      // constraints
      this->augmented_system.rhs[this->number_variables + j] = -this->barrier_constraints[j];
   }
   DEBUG << "RHS: "; print_vector(DEBUG, this->augmented_system.rhs, 0, this->number_variables + this->number_constraints); DEBUG << "\n";
}

void InteriorPoint::compute_lower_bound_dual_direction(const std::vector<double>& solution) {
   clear(this->lower_delta_z);
   for (size_t i: this->lower_bounded_variables) {
      const double distance_to_bound = this->primal_iterate[i] - this->variables_bounds[i].lb;
      this->lower_delta_z[i] = (this->barrier_parameter - solution[i] * this->lower_bound_multipliers[i]) / distance_to_bound -
                               this->lower_bound_multipliers[i];
   }
}

void InteriorPoint::compute_upper_bound_dual_direction(const std::vector<double>& solution) {
   clear(this->upper_delta_z);
   for (size_t i: this->upper_bounded_variables) {
      const double distance_to_bound = this->primal_iterate[i] - this->variables_bounds[i].ub;
      this->upper_delta_z[i] = (this->barrier_parameter - solution[i] * this->upper_bound_multipliers[i]) / distance_to_bound -
                               this->upper_bound_multipliers[i];
   }
}

void InteriorPoint::generate_direction(const Problem& problem, const Iterate& current_iterate) {
   // retrieve +Δλ (Nocedal p590)
   for (size_t j = this->number_variables; j < this->augmented_system.solution.size(); j++) {
      this->augmented_system.solution[j] = -this->augmented_system.solution[j];
   }

   // "fraction to boundary" rule for primal variables and constraints multipliers
   const double tau = std::max(this->parameters.tau_min, 1. - this->barrier_parameter);
   const double primal_step_length = this->primal_fraction_to_boundary(this->augmented_system.solution, tau);
   for (size_t i = 0; i < this->number_variables; i++) {
      this->direction.x[i] = primal_step_length * this->augmented_system.solution[i];
   }
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->direction.multipliers.constraints[j] = primal_step_length * this->augmented_system.solution[this->number_variables + j];
   }

   // compute bound multiplier direction Δz
   this->compute_lower_bound_dual_direction(this->augmented_system.solution);
   this->compute_upper_bound_dual_direction(this->augmented_system.solution);

   // "fraction to boundary" rule for bound multipliers
   const double dual_step_length = this->dual_fraction_to_boundary(tau);
   for (size_t i = 0; i < this->number_variables; i++) {
      this->direction.multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_step_length * this->lower_delta_z[i];
      this->direction.multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_step_length * this->upper_delta_z[i];
   }

   this->direction.norm = norm_inf(direction.x, 0, this->number_variables);
   // evaluate the barrier objective
   this->direction.objective = this->compute_barrier_directional_derivative(direction.x);
   this->print_solution(problem, primal_step_length, dual_step_length);
}

double InteriorPoint::compute_KKT_error_scaling(const Iterate& current_iterate) const {
   // KKT error
   const double norm_1_constraint_multipliers = norm_1(current_iterate.multipliers.constraints);
   const double norm_1_bound_multipliers = norm_1(current_iterate.multipliers.lower_bounds) + norm_1(current_iterate.multipliers.upper_bounds);
   const double norm_1_multipliers = norm_1_constraint_multipliers + norm_1_bound_multipliers;
   const size_t total_size = this->number_variables + current_iterate.multipliers.constraints.size();
   const double sd = std::max(this->parameters.smax, norm_1_multipliers / static_cast<double>(total_size)) / this->parameters.smax;
   return sd;
}

double InteriorPoint::compute_central_complementarity_error(const Iterate& iterate) const {
   // variable bound constraints
   const auto residual_function = [&](size_t i) {
      double result = 0.;
      if (is_finite_lower_bound(this->variables_bounds[i].lb)) {
         result += iterate.multipliers.lower_bounds[i] * (iterate.x[i] - this->variables_bounds[i].lb) - this->barrier_parameter;
      }
      if (is_finite_upper_bound(this->variables_bounds[i].ub)) {
         result += iterate.multipliers.upper_bounds[i] * (iterate.x[i] - this->variables_bounds[i].ub) - this->barrier_parameter;
      }
      return result;
   };

   // scaling
   const double bound_multipliers_norm = norm_1(iterate.multipliers.lower_bounds) + norm_1(iterate.multipliers.upper_bounds);
   const double sc = std::max(this->parameters.smax, bound_multipliers_norm / static_cast<double>(this->number_variables)) / this->parameters.smax;
   return norm_1(residual_function, this->number_variables) / sc;
}

void InteriorPoint::set_current_iterate(const Iterate& iterate) {
   copy_from(this->primal_iterate, iterate.x);
   copy_from(this->lower_bound_multipliers, iterate.multipliers.lower_bounds);
   copy_from(this->upper_bound_multipliers, iterate.multipliers.upper_bounds);
}

void InteriorPoint::register_accepted_iterate(Iterate& iterate) {
   // rescale the bound multipliers (Eq (16) in Ipopt paper)
   for (size_t i: this->lower_bounded_variables) {
      const double coefficient = this->barrier_parameter / (iterate.x[i] - this->variables_bounds[i].lb);
      const double lb = coefficient / this->parameters.k_sigma;
      const double ub = coefficient * this->parameters.k_sigma;
      assert(lb <= ub && "IPM bound multiplier reset: the bounds are in the wrong order");
      iterate.multipliers.lower_bounds[i] = std::max(std::min(iterate.multipliers.lower_bounds[i], ub), lb);
   }
   for (size_t i: this->upper_bounded_variables) {
      const double coefficient = this->barrier_parameter / (iterate.x[i] - this->variables_bounds[i].ub);
      const double lb = coefficient * this->parameters.k_sigma;
      const double ub = coefficient / this->parameters.k_sigma;
      assert(lb <= ub && "IPM bound multiplier reset: the bounds are in the wrong order");
      iterate.multipliers.upper_bounds[i] = std::max(std::min(iterate.multipliers.upper_bounds[i], ub), lb);
   }
}

size_t InteriorPoint::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}

void InteriorPoint::print_solution(const Problem& problem, double primal_step_length, double dual_step_length) const {
   DEBUG << "IPM solution:\n";
   DEBUG << "Δx: "; print_vector(DEBUG, this->augmented_system.solution, 0, problem.number_variables);
   DEBUG << "Δs: "; print_vector(DEBUG, this->augmented_system.solution, problem.number_variables, problem.inequality_constraints.size());
   if (this->number_variables > problem.number_variables + problem.inequality_constraints.size()) {
      DEBUG << "Δe: "; print_vector(DEBUG, this->augmented_system.solution, problem.number_variables + problem.inequality_constraints.size(),
            this->number_variables - (problem.number_variables + problem.inequality_constraints.size()));
   }
   DEBUG << "Δλ: "; print_vector(DEBUG, this->augmented_system.solution, this->number_variables, problem.number_constraints);
   DEBUG << "Δz_L: "; print_vector(DEBUG, this->lower_delta_z, 0, this->number_variables);
   DEBUG << "Δz_U: "; print_vector(DEBUG, this->upper_delta_z, 0, this->number_variables);
   DEBUG << "primal length = " << primal_step_length << "\n";
   DEBUG << "dual length = " << dual_step_length << "\n\n";
}