#include <cmath>
#include "BarrierSubproblem.hpp"
#include "solvers/linear/LinearSolverFactory.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"
#include "optimization/Preprocessing.hpp"

BarrierSubproblem::BarrierSubproblem(const Problem& problem, size_t max_number_variables, const Options& options):
      Subproblem(max_number_variables, // max_number_variables
            problem.number_constraints, SOC_UPON_REJECTION, true, norm_from_string(options.at("residual_norm"))),
      augmented_system(options.at("sparse_format"), max_number_variables + problem.number_constraints,
            problem.get_hessian_maximum_number_nonzeros()
            + max_number_variables + problem.number_constraints /* regularization */
            + 2 * max_number_variables /* diagonal barrier terms */
            + max_number_variables * problem.number_constraints /* Jacobian */,
            stod(options.at("LS_regularization_failure_threshold"))),
      barrier_parameter(std::stod(options.at("initial_barrier_parameter"))),
      previous_barrier_parameter(std::stod(options.at("initial_barrier_parameter"))),
      tolerance(std::stod(options.at("tolerance"))),
      // if no trust region is used, the problem should be convexified. However, the inertia of the augmented matrix will be corrected later
      hessian_model(HessianModelFactory::create(options.at("hessian_model"), max_number_variables, problem.get_hessian_maximum_number_nonzeros(),
            false, options)),
      linear_solver(LinearSolverFactory::create(options.at("linear_solver"), max_number_variables + problem.number_constraints,
            problem.get_hessian_maximum_number_nonzeros()
            + max_number_variables + problem.number_constraints /* regularization */
            + 2 * max_number_variables /* diagonal barrier terms */
            + max_number_variables * problem.number_constraints /* Jacobian */)),
      parameters({stod(options.at("tau_min")),
            stod(options.at("k_sigma")),
            stod(options.at("smax")),
            stod(options.at("k_mu")),
            stod(options.at("theta_mu")),
            stod(options.at("k_epsilon")),
            stod(options.at("barrier_update_fraction")),
            stod(options.at("regularization_barrier_exponent"))}),
      default_multiplier(std::stod(options.at("default_multiplier"))),
      lower_delta_z(max_number_variables), upper_delta_z(max_number_variables) {
   assert(problem.inequality_constraints.empty() && "The problem has inequality constraints. Create an instance of SlackReformulation");
   // register the variables bounds
   for (size_t i = 0; i < problem.number_variables; i++) {
      this->current_variable_bounds[i] = {problem.get_variable_lower_bound(i), problem.get_variable_upper_bound(i)};
   }
}

inline void BarrierSubproblem::initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) {
   statistics.add_column("barrier param.", Statistics::double_width, 8);

   // make the initial point strictly feasible wrt the bounds
   for (size_t i = 0; i < problem.number_variables; i++) {
      const Range bounds = {problem.get_variable_lower_bound(i), problem.get_variable_upper_bound(i)};
      first_iterate.x[i] = Subproblem::push_variable_to_interior(first_iterate.x[i], bounds);
   }

   // set the bound multipliers
   for (size_t i: problem.lower_bounded_variables) {
      first_iterate.multipliers.lower_bounds[i] = this->default_multiplier;
   }
   for (size_t i: problem.upper_bounded_variables) {
      first_iterate.multipliers.upper_bounds[i] = -this->default_multiplier;
   }

   // compute least-square multipliers
   if (problem.is_constrained()) {
      this->augmented_system.matrix->dimension = problem.number_variables + problem.number_constraints;
      Preprocessing::compute_least_square_multipliers(problem, *this->augmented_system.matrix, this->augmented_system.rhs, *this->linear_solver,
            first_iterate, first_iterate.multipliers.constraints);
   }
}

void BarrierSubproblem::build_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier,
      double trust_region_radius) {
   // update the barrier parameter if the current iterate solves the subproblem
   this->update_barrier_parameter(problem, current_iterate);

   // constraint Jacobian
   current_iterate.evaluate_constraint_jacobian(problem);

   // build a model of the objective scaled by the objective multiplier
   this->build_objective_model(problem, current_iterate, objective_multiplier);

   // variables and bounds
   this->set_current_variable_bounds(problem, current_iterate, trust_region_radius);
}

void BarrierSubproblem::build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // if we're building the feasibility subproblem, temporarily update the objective multiplier
   if (objective_multiplier == 0.) {
      this->solving_feasibility_problem = true;
      this->previous_barrier_parameter = this->barrier_parameter;
      this->barrier_parameter = std::max(this->barrier_parameter, norm_inf(current_iterate.constraints));
      DEBUG << "Barrier parameter mu temporarily updated to " << this->barrier_parameter << "\n";
      this->subproblem_definition_changed = true;
   }
   else {
      this->solving_feasibility_problem = false;
   }

   // evaluate the Hessian
   this->hessian_model->evaluate(problem, current_iterate.x, objective_multiplier, current_iterate.multipliers.constraints);

   // objective gradient
   this->set_scaled_objective_gradient(problem, current_iterate, objective_multiplier);
   for (size_t i: problem.lower_bounded_variables) {
      this->objective_gradient.insert(i, -this->barrier_parameter / (current_iterate.x[i] - this->current_variable_bounds[i].lb));
   }
   for (size_t i: problem.upper_bounded_variables) {
      this->objective_gradient.insert(i, -this->barrier_parameter / (current_iterate.x[i] - this->current_variable_bounds[i].ub));
   }
   this->augmented_system.matrix->reset();
}

Direction BarrierSubproblem::solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   // set up the augmented system (with the correct inertia)
   this->hessian_model->adjust_number_variables(problem.number_variables);
   this->assemble_augmented_system(problem, current_iterate);

   // compute the solution (Δx, -Δλ)
   this->augmented_system.solve(*this->linear_solver);
   assert(this->direction.status == OPTIMAL && "The barrier subproblem was not solved to optimality");
   this->number_subproblems_solved++;

   // generate IPM direction
   this->generate_direction(problem, current_iterate);
   statistics.add_statistic("barrier param.", this->barrier_parameter);

   // determine if the direction is a "small direction" (Section 3.9 of the Ipopt paper) TODO
   bool is_small_direction = this->is_small_direction(problem, current_iterate, this->direction);
   if (is_small_direction) {
      DEBUG << "This is a small direction\n";
   }
   return this->direction;
}

void BarrierSubproblem::assemble_augmented_system(const Problem& problem, const Iterate& current_iterate) {
   // assemble, factorize and regularize the KKT matrix
   this->assemble_augmented_matrix(problem, current_iterate);
   this->augmented_system.factorize_matrix(problem, *this->linear_solver);
   this->augmented_system.regularize_matrix(problem, *this->linear_solver, problem.number_variables, problem.number_constraints,
         std::pow(this->barrier_parameter, this->parameters.regularization_barrier_exponent));
   auto[number_pos, number_neg, number_zero] = this->linear_solver->get_inertia();
   assert(number_pos == problem.number_variables && number_neg == problem.number_constraints && number_zero == 0);

   // right-hand side
   this->generate_augmented_rhs(problem, current_iterate);
}

Direction BarrierSubproblem::compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) {
   DEBUG << "\nEntered SOC computation\n";
   // modify the RHS by adding the values of the constraints
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->augmented_system.rhs[problem.number_variables + j] -= trial_iterate.constraints[j];
   }
   DEBUG << "SOC RHS: "; print_vector(DEBUG, this->augmented_system.rhs, 0, problem.number_variables + problem.number_constraints);

   // compute the solution (Δx, -Δλ)
   this->augmented_system.solve(*this->linear_solver);
   this->number_subproblems_solved++;

   // generate IPM direction
   this->generate_direction(problem, trial_iterate);
   return this->direction;
}

double BarrierSubproblem::get_proximal_coefficient() const {
   return std::sqrt(this->barrier_parameter)/2.;
}

PredictedReductionModel BarrierSubproblem::generate_predicted_reduction_model(const Problem& /*problem*/, const Iterate& /*current_iterate*/,
      const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() {
      return [=](double step_length) {
         return -step_length * direction.objective;
      };
   });
}

double BarrierSubproblem::compute_optimality_measure(const Problem& problem, Iterate& iterate) {
   return this->evaluate_barrier_function(problem, iterate);
}

void BarrierSubproblem::update_barrier_parameter(const Problem& problem, const Iterate& current_iterate) {
   // scaled error terms
   const double sd = this->compute_KKT_error_scaling(problem, current_iterate);
   const double KKTerror = current_iterate.nonlinear_errors.stationarity / sd;
   const double central_complementarity_error = this->compute_central_complementarity_error(problem, current_iterate);
   const double error = std::max({KKTerror, current_iterate.nonlinear_errors.constraints, central_complementarity_error});
   DEBUG << "KKT error for barrier subproblem is " << error << "\n";

   // update of the barrier parameter (Eq. 7 in Ipopt paper)
   const double tolerance_fraction = this->tolerance / this->parameters.barrier_update_fraction;
   while (error <= this->parameters.k_epsilon * this->barrier_parameter && tolerance_fraction < this->barrier_parameter) {
      this->barrier_parameter = std::max(tolerance_fraction, std::min(this->parameters.k_mu * this->barrier_parameter,
            std::pow(this->barrier_parameter, this->parameters.theta_mu)));
      DEBUG << "Barrier parameter mu updated to " << this->barrier_parameter << "\n";
      // signal the redefinition of the problem to the globalization strategy
      this->subproblem_definition_changed = true;
   }
}

bool BarrierSubproblem::is_small_direction(const Problem& problem, const Iterate& current_iterate, const Direction& direction) {
   const auto relative_measure_function = [&](size_t i) {
      return direction.x[i]/(1 + current_iterate.x[i]);
   };
   const double machine_epsilon = std::numeric_limits<double>::epsilon();
   return (norm_inf(relative_measure_function, problem.number_variables) < 10. * machine_epsilon);
}

void BarrierSubproblem::set_current_variable_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   // apply the trust region
   for (size_t i = 0; i < problem.number_variables; i++) {
      double lb = std::max(current_iterate.x[i] - trust_region_radius, problem.get_variable_lower_bound(i));
      double ub = std::min(current_iterate.x[i] + trust_region_radius, problem.get_variable_upper_bound(i));
      this->current_variable_bounds[i] = {lb, ub};
   }
}

double BarrierSubproblem::compute_barrier_directional_derivative(const std::vector<double>& solution) {
   return dot(solution, this->objective_gradient);
}

double BarrierSubproblem::evaluate_barrier_function(const Problem& problem, Iterate& iterate) {
   double objective = 0.;
   // bound constraints
   for (size_t i: problem.lower_bounded_variables) {
      objective -= std::log(iterate.x[i] - this->current_variable_bounds[i].lb);
   }
   for (size_t i: problem.upper_bounded_variables) {
      objective -= std::log(this->current_variable_bounds[i].ub - iterate.x[i]);
   }
   objective *= this->barrier_parameter;
   if (!this->solving_feasibility_problem) {
      // original objective
      iterate.evaluate_objective(problem);
      objective += iterate.objective;
   }
   return objective;
}

double BarrierSubproblem::primal_fraction_to_boundary(const Problem& problem, const Iterate& current_iterate, const std::vector<double>& ipm_solution,
      double tau) {
   double primal_length = 1.;
   for (size_t i: problem.lower_bounded_variables) {
      if (ipm_solution[i] < 0.) {
         double trial_alpha_xi = -tau * (current_iterate.x[i] - this->current_variable_bounds[i].lb) / ipm_solution[i];
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   for (size_t i: problem.upper_bounded_variables) {
      if (0. < ipm_solution[i]) {
         double trial_alpha_xi = -tau * (current_iterate.x[i] - this->current_variable_bounds[i].ub) / ipm_solution[i];
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   return primal_length;
}

double BarrierSubproblem::dual_fraction_to_boundary(const Problem& problem, const Iterate& current_iterate, double tau) {
   double dual_length = 1.;
   for (size_t i = 0; i < problem.number_variables; i++) {
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

void BarrierSubproblem::assemble_augmented_matrix(const Problem& problem, const Iterate& current_iterate) {
   this->augmented_system.matrix->dimension = problem.number_variables + problem.number_constraints;
   // copy the Lagrangian Hessian in the top left block
   size_t current_column = 0;
   this->hessian_model->hessian->for_each([&](size_t i, size_t j, double entry) {
      for (size_t column = current_column; column < j; column++) {
         this->augmented_system.matrix->finalize(column);
         current_column++;
      }
      this->augmented_system.matrix->insert(entry, i, j);
   });

   // diagonal barrier terms
   for (size_t i: problem.lower_bounded_variables) {
      const double diagonal_term = current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - this->current_variable_bounds[i].lb);
      this->augmented_system.matrix->insert(diagonal_term, i, i);
   }
   for (size_t i: problem.upper_bounded_variables) {
      const double diagonal_term = current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - this->current_variable_bounds[i].ub);
      this->augmented_system.matrix->insert(diagonal_term, i, i);
   }

   // Jacobian of general constraints
   for (size_t j = 0; j < problem.number_constraints; j++) {
      current_iterate.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         this->augmented_system.matrix->insert(derivative, i, problem.number_variables + j);
      });
      this->augmented_system.matrix->finalize(j);
   }
}

void BarrierSubproblem::generate_augmented_rhs(const Problem& problem, const Iterate& current_iterate) {
   // generate the right-hand side
   initialize_vector(this->augmented_system.rhs, 0.);

   // barrier objective gradient
   this->objective_gradient.for_each([&](size_t i, double derivative) {
      this->augmented_system.rhs[i] = -derivative;
   });

   // constraint: evaluations and gradients
   for (size_t j = 0; j < problem.number_constraints; j++) {
      // Lagrangian
      if (current_iterate.multipliers.constraints[j] != 0.) {
         current_iterate.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
            this->augmented_system.rhs[i] += current_iterate.multipliers.constraints[j] * derivative;
         });
      }
      // constraints
      this->augmented_system.rhs[problem.number_variables + j] = -current_iterate.constraints[j];
   }
   DEBUG << "RHS: "; print_vector(DEBUG, this->augmented_system.rhs, 0, problem.number_variables + problem.number_constraints); DEBUG << "\n";
}

void BarrierSubproblem::compute_lower_bound_dual_direction(const Problem& problem, const Iterate& current_iterate, const std::vector<double>& solution) {
   initialize_vector(this->lower_delta_z, 0.);
   for (size_t i: problem.lower_bounded_variables) {
      const double distance_to_bound = current_iterate.x[i] - this->current_variable_bounds[i].lb;
      this->lower_delta_z[i] = (this->barrier_parameter - solution[i] * current_iterate.multipliers.lower_bounds[i]) / distance_to_bound -
            current_iterate.multipliers.lower_bounds[i];
   }
}

void BarrierSubproblem::compute_upper_bound_dual_direction(const Problem& problem, const Iterate& current_iterate, const std::vector<double>& solution) {
   initialize_vector(this->upper_delta_z, 0.);
   for (size_t i: problem.upper_bounded_variables) {
      const double distance_to_bound = current_iterate.x[i] - this->current_variable_bounds[i].ub;
      this->upper_delta_z[i] = (this->barrier_parameter - solution[i] * current_iterate.multipliers.upper_bounds[i]) / distance_to_bound -
            current_iterate.multipliers.upper_bounds[i];
   }
}

void BarrierSubproblem::generate_direction(const Problem& problem, const Iterate& current_iterate) {
   // retrieve +Δλ (Nocedal p590)
   for (size_t j = problem.number_variables; j < this->augmented_system.solution.size(); j++) {
      this->augmented_system.solution[j] = -this->augmented_system.solution[j];
   }

   // "fraction to boundary" rule for primal variables and constraints multipliers
   const double tau = std::max(this->parameters.tau_min, 1. - this->barrier_parameter);
   const double primal_step_length = this->primal_fraction_to_boundary(problem, current_iterate, this->augmented_system.solution, tau);
   for (size_t i = 0; i < problem.number_variables; i++) {
      this->direction.x[i] = primal_step_length * this->augmented_system.solution[i];
   }
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->direction.multipliers.constraints[j] = primal_step_length * this->augmented_system.solution[problem.number_variables + j];
   }

   // compute bound multiplier direction Δz
   this->compute_lower_bound_dual_direction(problem, current_iterate, this->augmented_system.solution);
   this->compute_upper_bound_dual_direction(problem, current_iterate, this->augmented_system.solution);

   // "fraction to boundary" rule for bound multipliers
   const double dual_step_length = this->dual_fraction_to_boundary(problem, current_iterate, tau);
   for (size_t i = 0; i < problem.number_variables; i++) {
      this->direction.multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_step_length * this->lower_delta_z[i];
      this->direction.multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_step_length * this->upper_delta_z[i];
   }

   this->direction.norm = norm_inf(direction.x, 0, problem.number_variables);
   // evaluate the barrier objective
   this->direction.objective = this->compute_barrier_directional_derivative(direction.x);
   this->print_solution(problem, primal_step_length, dual_step_length);
}

double BarrierSubproblem::compute_KKT_error_scaling(const Problem& problem, const Iterate& current_iterate) const {
   // KKT error
   const double norm_1_constraint_multipliers = norm_1(current_iterate.multipliers.constraints);
   const double norm_1_bound_multipliers = norm_1(current_iterate.multipliers.lower_bounds) + norm_1(current_iterate.multipliers.upper_bounds);
   const double norm_1_multipliers = norm_1_constraint_multipliers + norm_1_bound_multipliers;
   const size_t total_size = problem.number_variables + problem.number_constraints;
   const double sd = std::max(this->parameters.smax, norm_1_multipliers / static_cast<double>(total_size)) / this->parameters.smax;
   return sd;
}

double BarrierSubproblem::compute_central_complementarity_error(const Problem& problem, const Iterate& iterate) const {
   // variable bounds TODO use problem.lower_bounded_variables
   const auto residual_function = [&](size_t i) {
      double result = 0.;
      if (is_finite(this->current_variable_bounds[i].lb)) {
         result += iterate.multipliers.lower_bounds[i] * (iterate.x[i] - this->current_variable_bounds[i].lb) - this->barrier_parameter;
      }
      if (is_finite(this->current_variable_bounds[i].ub)) {
         result += iterate.multipliers.upper_bounds[i] * (iterate.x[i] - this->current_variable_bounds[i].ub) - this->barrier_parameter;
      }
      return result;
   };

   // scaling
   const double bound_multipliers_norm = norm_1(iterate.multipliers.lower_bounds) + norm_1(iterate.multipliers.upper_bounds);
   const double sc = std::max(this->parameters.smax, bound_multipliers_norm / static_cast<double>(problem.number_variables)) / this->parameters.smax;
   return norm_1(residual_function, problem.number_variables) / sc;
}

void BarrierSubproblem::register_accepted_iterate(const Problem& problem, Iterate& iterate) {
   if (this->solving_feasibility_problem) {
       this->barrier_parameter = this->previous_barrier_parameter;
       this->solving_feasibility_problem = false;
   }
   if (this->solving_feasibility_problem) {
      // compute least-square multipliers TODO
   }

   // rescale the bound multipliers (Eq. 16 in Ipopt paper)
   for (size_t i: problem.lower_bounded_variables) {
      const double coefficient = this->barrier_parameter / (iterate.x[i] - this->current_variable_bounds[i].lb);
      const double lb = coefficient / this->parameters.k_sigma;
      const double ub = coefficient * this->parameters.k_sigma;
      assert(lb <= ub && "IPM lower bound multiplier reset: the bounds are in the wrong order");
      iterate.multipliers.lower_bounds[i] = std::max(std::min(iterate.multipliers.lower_bounds[i], ub), lb);
   }
   for (size_t i: problem.upper_bounded_variables) {
      const double coefficient = this->barrier_parameter / (iterate.x[i] - this->current_variable_bounds[i].ub);
      const double lb = coefficient * this->parameters.k_sigma;
      const double ub = coefficient / this->parameters.k_sigma;
      assert(lb <= ub && "IPM upper bound multiplier reset: the bounds are in the wrong order");
      iterate.multipliers.upper_bounds[i] = std::max(std::min(iterate.multipliers.upper_bounds[i], ub), lb);
   }
}

size_t BarrierSubproblem::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}

void BarrierSubproblem::print_solution(const Problem& problem, double primal_step_length, double dual_step_length) const {
   DEBUG << "IPM solution:\n";
   DEBUG << "Δx: "; print_vector(DEBUG, this->augmented_system.solution, 0, problem.number_variables);
   if (problem.get_number_original_variables() < problem.number_variables) {
      DEBUG << "Δe: "; print_vector(DEBUG, this->augmented_system.solution, problem.number_variables, problem.number_variables - problem.number_variables);
   }
   DEBUG << "Δλ: "; print_vector(DEBUG, this->augmented_system.solution, problem.number_variables, problem.number_constraints);
   DEBUG << "Δz_L: "; print_vector(DEBUG, this->lower_delta_z, 0, problem.number_variables);
   DEBUG << "Δz_U: "; print_vector(DEBUG, this->upper_delta_z, 0, problem.number_variables);
   DEBUG << "primal length = " << primal_step_length << "\n";
   DEBUG << "dual length = " << dual_step_length << "\n";
}

void BarrierSubproblem::set_initial_point(const std::vector<double>& /*initial_point*/) {
   // do nothing
}