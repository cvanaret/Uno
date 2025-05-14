// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "PrimalDualInteriorPointMethod.hpp"
#include "PrimalDualInteriorPointProblem.hpp"
#include "ingredients/constraint_relaxation_strategies/l1RelaxedProblem.hpp"
#include "ingredients/regularization_strategies/Inertia.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "layers/SubproblemLayer.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "preprocessing/Preprocessing.hpp"
#include "symbolic/Expression.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   PrimalDualInteriorPointMethod::PrimalDualInteriorPointMethod(const Options& options):
         InequalityHandlingMethod(),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"))),
         barrier_parameter_update_strategy(options),
         previous_barrier_parameter(options.get_double("barrier_initial_parameter")),
         default_multiplier(options.get_double("barrier_default_multiplier")),
         parameters({
               options.get_double("barrier_tau_min"),
               options.get_double("barrier_k_sigma"),
               options.get_double("barrier_regularization_exponent"),
               options.get_double("barrier_small_direction_factor"),
               options.get_double("barrier_push_variable_to_interior_k1"),
               options.get_double("barrier_push_variable_to_interior_k2")
         }),
         least_square_multiplier_max_norm(options.get_double("least_square_multiplier_max_norm")),
         damping_factor(options.get_double("barrier_damping_factor")),
         l1_constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")) {
   }

   void PrimalDualInteriorPointMethod::initialize(const OptimizationProblem& problem, const HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy) {
      const PrimalDualInteriorPointProblem barrier_problem(problem, this->barrier_parameter());
      this->objective_gradient.reserve(barrier_problem.number_objective_gradient_nonzeros());
      this->constraints.resize(barrier_problem.number_constraints);
      this->constraint_jacobian.resize(barrier_problem.number_constraints, barrier_problem.number_variables);
      const size_t number_hessian_nonzeros = barrier_problem.number_hessian_nonzeros(hessian_model);
      const size_t number_augmented_system_nonzeros = number_hessian_nonzeros + barrier_problem.number_jacobian_nonzeros() +
         barrier_problem.number_variables + barrier_problem.number_constraints; // primal-dual regularization,
      this->hessian = SymmetricMatrix<size_t, double>(barrier_problem.number_variables, number_hessian_nonzeros, false, "COO");
      this->linear_solver->initialize_memory(barrier_problem.number_variables + barrier_problem.number_constraints,
         number_augmented_system_nonzeros);
      regularization_strategy.initialize_memory(barrier_problem, hessian_model);
      this->augmented_system.matrix = SymmetricMatrix<size_t, double>(barrier_problem.number_variables + barrier_problem.number_constraints,
         number_augmented_system_nonzeros,
         true, "COO");
      this->augmented_system.rhs.resize(barrier_problem.number_variables + barrier_problem.number_constraints);
      this->augmented_system.solution.resize(barrier_problem.number_variables + barrier_problem.number_constraints);
   }

   void PrimalDualInteriorPointMethod::initialize_statistics(Statistics& statistics, const Options& options) {
      statistics.add_column("regulariz", Statistics::double_width - 4, options.get_int("statistics_regularization_column_order"));
      statistics.add_column("barrier", Statistics::double_width - 5, options.get_int("statistics_barrier_parameter_column_order"));
   }

   void PrimalDualInteriorPointMethod::generate_initial_iterate(const OptimizationProblem& first_reformulation, Iterate& initial_iterate) {
      // TODO: enforce linear constraints at initial point
      //if (options.get_bool("enforce_linear_constraints")) {
      //   Preprocessing::enforce_linear_constraints(problem.model, initial_iterate.primals, initial_iterate.multipliers, this->solver);
      //}

      const PrimalDualInteriorPointProblem barrier_problem(first_reformulation, this->barrier_parameter());

      // add the slacks to the initial iterate
      initial_iterate.set_number_variables(barrier_problem.number_variables);
      // make the initial point strictly feasible wrt the bounds
      for (size_t variable_index: Range(first_reformulation.number_variables)) {
         initial_iterate.primals[variable_index] = PrimalDualInteriorPointMethod::push_variable_to_interior(initial_iterate.primals[variable_index],
            first_reformulation.variable_lower_bound(variable_index), first_reformulation.variable_upper_bound(variable_index));
      }

      // set the slack variables (if any)
      if (!first_reformulation.get_inequality_constraints().empty()) {
         // evaluate the constraints at the original point
         initial_iterate.evaluate_constraints(first_reformulation.model);

         // set the slacks to the constraint values
         size_t slack_index = first_reformulation.number_variables;
         for (const auto inequality_constraint_index: first_reformulation.get_inequality_constraints()) {
            initial_iterate.primals[slack_index] = PrimalDualInteriorPointMethod::push_variable_to_interior(
               initial_iterate.evaluations.constraints[inequality_constraint_index],
               first_reformulation.variable_lower_bound(slack_index), first_reformulation.variable_upper_bound(slack_index));
            slack_index++;
         }
         // since the slacks have been set, the function evaluations should also be updated
         initial_iterate.is_objective_gradient_computed = false;
         initial_iterate.are_constraints_computed = false;
         initial_iterate.is_constraint_jacobian_computed = false;
      }

      // set the bound multipliers
      for (const size_t variable_index: first_reformulation.get_lower_bounded_variables()) {
         initial_iterate.multipliers.lower_bounds[variable_index] = this->default_multiplier;
      }
      for (const size_t variable_index: first_reformulation.get_upper_bounded_variables()) {
         initial_iterate.multipliers.upper_bounds[variable_index] = -this->default_multiplier;
      }
      size_t slack_index = first_reformulation.number_variables;
      for (size_t inequality_constraint_index: first_reformulation.get_inequality_constraints()) {
         if (is_finite(first_reformulation.constraint_lower_bound(inequality_constraint_index))) {
            initial_iterate.multipliers.lower_bounds[slack_index] = this->default_multiplier;
         }
         if (is_finite(first_reformulation.constraint_upper_bound(inequality_constraint_index))) {
            initial_iterate.multipliers.lower_bounds[slack_index] = -this->default_multiplier;
         }
         slack_index++;
      }

      // compute least-square multipliers
      if (0 < first_reformulation.number_constraints) {
         //this->compute_least_square_multipliers(first_reformulation, initial_iterate, initial_iterate.multipliers.constraints);
      }
   }

   double PrimalDualInteriorPointMethod::barrier_parameter() const {
      return this->barrier_parameter_update_strategy.get_barrier_parameter();
   }

   double PrimalDualInteriorPointMethod::push_variable_to_interior(double variable_value, double lower_bound, double upper_bound) const {
      const double range = upper_bound - lower_bound;
      const double perturbation_lb = std::min(this->parameters.push_variable_to_interior_k1 * std::max(1., std::abs(lower_bound)),
         this->parameters.push_variable_to_interior_k2 * range);
      const double perturbation_ub = std::min(this->parameters.push_variable_to_interior_k1 * std::max(1., std::abs(upper_bound)),
         this->parameters.push_variable_to_interior_k2 * range);
      variable_value = std::max(variable_value, lower_bound + perturbation_lb);
      variable_value = std::min(variable_value, upper_bound - perturbation_ub);
      return variable_value;
   }

   void PrimalDualInteriorPointMethod::evaluate_functions(Statistics& statistics, const PrimalDualInteriorPointProblem& barrier_problem,
         Iterate& current_iterate, const Multipliers& current_multipliers, SubproblemLayer& subproblem_layer, const WarmstartInformation& warmstart_information) {
      // barrier objective gradient
      if (warmstart_information.objective_changed) {
         barrier_problem.evaluate_objective_gradient(current_iterate, this->objective_gradient);
      }

      // constraints and Jacobian
      if (warmstart_information.constraints_changed) {
         barrier_problem.evaluate_constraints(current_iterate, this->constraints);
         barrier_problem.evaluate_constraint_jacobian(current_iterate, this->constraint_jacobian);
      }

      // barrier Lagrangian Hessian
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         barrier_problem.evaluate_lagrangian_hessian(statistics, *subproblem_layer.hessian_model, current_iterate.primals, current_multipliers, this->hessian);
      }
   }

   void PrimalDualInteriorPointMethod::solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, SubproblemLayer& subproblem_layer, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("The interior-point subproblem has a trust region. This is not implemented yet");
      }

      // possibly update the barrier parameter
      const auto& residuals = this->solving_feasibility_problem ? current_iterate.feasibility_residuals : current_iterate.residuals;
      if (!this->first_feasibility_iteration) {
         this->update_barrier_parameter(problem, current_iterate, current_multipliers, residuals);
      }
      else {
         this->first_feasibility_iteration = false;
      }
      statistics.set("barrier", this->barrier_parameter());

      // create a barrier problem
      const PrimalDualInteriorPointProblem barrier_problem(problem, this->barrier_parameter());
      direction.set_dimensions(barrier_problem.number_variables, barrier_problem.number_constraints);

      // evaluate the functions at the current iterate
      this->evaluate_functions(statistics, barrier_problem, current_iterate, current_multipliers, subproblem_layer, warmstart_information);

      // compute the primal-dual solution
      this->assemble_augmented_system(statistics, problem, current_multipliers, subproblem_layer);
      this->augmented_system.solve(*this->linear_solver);
      assert(direction.status == SubproblemStatus::OPTIMAL && "The primal-dual perturbed subproblem was not solved to optimality");
      this->number_subproblems_solved++;

      this->assemble_primal_dual_direction(barrier_problem, current_iterate.primals, current_multipliers, direction.primals, direction.multipliers);
      direction.subproblem_objective = this->evaluate_subproblem_objective(direction);
   }

   double PrimalDualInteriorPointMethod::hessian_quadratic_product(const Vector<double>& /*vector*/) const {
      return 0.; // TODO
   }

   void PrimalDualInteriorPointMethod::assemble_augmented_system(Statistics& statistics, const OptimizationProblem& problem,
         const Multipliers& current_multipliers, SubproblemLayer& subproblem_layer) {
      // assemble and factorize the augmented matrix
      this->augmented_system.assemble_matrix(this->hessian, this->constraint_jacobian, problem.number_variables, problem.number_constraints);

      // regularize the augmented matrix
      const double dual_regularization_parameter = std::pow(this->barrier_parameter(), this->parameters.regularization_exponent);
      const Inertia expected_inertia{problem.number_variables, problem.number_constraints, 0};
      subproblem_layer.regularization_strategy->regularize_augmented_matrix(statistics, this->augmented_system.matrix,
         dual_regularization_parameter, expected_inertia, *this->linear_solver);

      // assemble the RHS
      this->assemble_augmented_rhs(current_multipliers, problem.number_variables, problem.number_constraints);
   }

   void PrimalDualInteriorPointMethod::initialize_feasibility_problem(const l1RelaxedProblem& /*problem*/, Iterate& current_iterate) {
      this->solving_feasibility_problem = true;
      this->first_feasibility_iteration = true;
      this->subproblem_definition_changed = true;

      // temporarily update the objective multiplier
      this->previous_barrier_parameter = this->barrier_parameter();
      const double new_barrier_parameter = std::max(this->barrier_parameter(), current_iterate.primal_feasibility);
      this->barrier_parameter_update_strategy.set_barrier_parameter(new_barrier_parameter);
      DEBUG << "Barrier parameter mu temporarily updated to " << this->barrier_parameter() << '\n';

      // set the bound multipliers
      /*
      for (const size_t variable_index: problem.get_lower_bounded_variables()) {
         current_iterate.feasibility_multipliers.lower_bounds[variable_index] = std::min(this->default_multiplier, problem.constraint_violation_coefficient);
      }
      for (const size_t variable_index: problem.get_upper_bounded_variables()) {
         current_iterate.feasibility_multipliers.upper_bounds[variable_index] = -this->default_multiplier;
      }
       */
   }

   // set the elastic variables of the current iterate
   void PrimalDualInteriorPointMethod::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) {
      DEBUG << "IPM: setting the elastic variables and their duals\n";

      // TODO set the bound multipliers
      for (const size_t variable_index: problem.get_lower_bounded_variables()) {
         current_iterate.feasibility_multipliers.lower_bounds[variable_index] = this->default_multiplier;
      }
      for (const size_t variable_index: problem.get_upper_bounded_variables()) {
         current_iterate.feasibility_multipliers.upper_bounds[variable_index] = -this->default_multiplier;
      }

      // c(x) - p + n = 0
      // analytical expression for p and n:
      // (mu_over_rho - jacobian_coefficient*this->barrier_constraints[j] + std::sqrt(radical))/2.
      // where jacobian_coefficient = -1 for p, +1 for n
      // Note: IPOPT uses a '+' sign because they define the Lagrangian as f(x) + \lambda^T c(x)
      const double mu = this->barrier_parameter();
      const auto elastic_setting_function = [&](Iterate& iterate, size_t constraint_index, size_t elastic_index, double jacobian_coefficient) {
         // precomputations
         const double constraint_j = this->constraints[constraint_index];
         const double rho = this->l1_constraint_violation_coefficient;
         const double mu_over_rho = mu / rho;
         const double radical = std::pow(constraint_j, 2) + std::pow(mu_over_rho, 2);
         const double sqrt_radical = std::sqrt(radical);

         iterate.primals[elastic_index] = (mu_over_rho - jacobian_coefficient * constraint_j + sqrt_radical) / 2.;
         iterate.feasibility_multipliers.lower_bounds[elastic_index] = mu / iterate.primals[elastic_index];
         iterate.feasibility_multipliers.upper_bounds[elastic_index] = 0.;
         assert(0. < iterate.primals[elastic_index] && "The elastic variable is not strictly positive.");
         assert(0. < iterate.feasibility_multipliers.lower_bounds[elastic_index] && "The elastic dual is not strictly positive.");
      };
      problem.set_elastic_variable_values(current_iterate, elastic_setting_function);
   }

   double PrimalDualInteriorPointMethod::proximal_coefficient() const {
      return std::sqrt(this->barrier_parameter());
   }

   void PrimalDualInteriorPointMethod::exit_feasibility_problem(const OptimizationProblem& /*problem*/, Iterate& /*trial_iterate*/) {
      //assert(this->solving_feasibility_problem && "The barrier subproblem did not know it was solving the feasibility problem.");
      this->barrier_parameter_update_strategy.set_barrier_parameter(this->previous_barrier_parameter);
      this->solving_feasibility_problem = false;
      // TODO computing the least-square multipliers messes up the factorization in the linear solver
      //this->compute_least_square_multipliers(problem, trial_iterate, trial_iterate.multipliers.constraints);
   }

   void PrimalDualInteriorPointMethod::set_auxiliary_measure(const Model& model, Iterate& iterate) {
      // auxiliary measure: barrier terms
      double barrier_terms = 0.;
      for (const size_t variable_index: model.get_lower_bounded_variables()) {
         barrier_terms -= std::log(iterate.primals[variable_index] - model.variable_lower_bound(variable_index));
      }
      for (const size_t variable_index: model.get_upper_bounded_variables()) {
         barrier_terms -= std::log(model.variable_upper_bound(variable_index) - iterate.primals[variable_index]);
      }
      // damping
      for (const size_t variable_index: model.get_single_lower_bounded_variables()) {
         barrier_terms += this->damping_factor*(iterate.primals[variable_index] - model.variable_lower_bound(variable_index));
      }
      for (const size_t variable_index: model.get_single_upper_bounded_variables()) {
         barrier_terms += this->damping_factor*(model.variable_upper_bound(variable_index) - iterate.primals[variable_index]);
      }
      barrier_terms *= this->barrier_parameter();
      assert(!std::isnan(barrier_terms) && "The auxiliary measure is not an number.");
      iterate.progress.auxiliary = barrier_terms;
   }

   double PrimalDualInteriorPointMethod::compute_predicted_auxiliary_reduction_model(const Model& model, const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const {
      const double directional_derivative = this->compute_barrier_term_directional_derivative(model, current_iterate, primal_direction);
      return step_length * (-directional_derivative);
      // }, "α*(μ*X^{-1} e^T d)"};
   }

   double PrimalDualInteriorPointMethod::compute_barrier_term_directional_derivative(const Model& model, const Iterate& current_iterate,
         const Vector<double>& primal_direction) const {
      double directional_derivative = 0.;
      for (const size_t variable_index: model.get_lower_bounded_variables()) {
         directional_derivative += -this->barrier_parameter() / (current_iterate.primals[variable_index] -
                                                                 model.variable_lower_bound(variable_index)) * primal_direction[variable_index];
      }
      for (const size_t variable_index: model.get_upper_bounded_variables()) {
         directional_derivative += -this->barrier_parameter() / (current_iterate.primals[variable_index] -
                                                                 model.variable_upper_bound(variable_index)) * primal_direction[variable_index];
      }
      // damping
      for (const size_t variable_index: model.get_single_lower_bounded_variables()) {
         directional_derivative += this->damping_factor * this->barrier_parameter() * primal_direction[variable_index];
      }
      for (const size_t variable_index: model.get_single_upper_bounded_variables()) {
         directional_derivative -= this->damping_factor * this->barrier_parameter() * primal_direction[variable_index];
      }
      return directional_derivative;
   }

   void PrimalDualInteriorPointMethod::update_barrier_parameter(const OptimizationProblem& problem, const Iterate& current_iterate,
         const Multipliers& current_multipliers, const DualResiduals& residuals) {
      const bool barrier_parameter_updated = this->barrier_parameter_update_strategy.update_barrier_parameter(problem, current_iterate,
            current_multipliers, residuals);
      // the barrier parameter may have been changed earlier when entering restoration
      this->subproblem_definition_changed = this->subproblem_definition_changed || barrier_parameter_updated;
   }

   // Section 3.9 in IPOPT paper
   bool PrimalDualInteriorPointMethod::is_small_step(const OptimizationProblem& problem, const Vector<double>& current_primals,
         const Vector<double>& direction_primals) const {
      const Range variables_range = Range(problem.number_variables);
      const VectorExpression relative_direction_size{variables_range, [&](size_t variable_index) {
         return direction_primals[variable_index] / (1 + std::abs(current_primals[variable_index]));
      }};
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      return (norm_inf(relative_direction_size) <= this->parameters.small_direction_factor * machine_epsilon);
   }

   double PrimalDualInteriorPointMethod::evaluate_subproblem_objective(const Direction& direction) const {
      const double linear_term = dot(direction.primals, this->objective_gradient);
      const double quadratic_term = this->hessian.quadratic_product(direction.primals, direction.primals) / 2.;
      return linear_term + quadratic_term;
   }

   // TODO use a single function for primal and dual fraction-to-boundary rules
   double PrimalDualInteriorPointMethod::primal_fraction_to_boundary(const OptimizationProblem& problem, const Vector<double>& current_primals,
         const Vector<double>& primal_direction, double tau) {
      double step_length = 1.;
      for (const size_t variable_index: problem.get_lower_bounded_variables()) {
         if (primal_direction[variable_index] < 0.) {
            double distance = -tau * (current_primals[variable_index] - problem.variable_lower_bound(variable_index)) / primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      for (const size_t variable_index: problem.get_upper_bounded_variables()) {
         if (0. < primal_direction[variable_index]) {
            double distance = -tau * (current_primals[variable_index] - problem.variable_upper_bound(variable_index)) / primal_direction[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      assert(0. < step_length && step_length <= 1. && "The primal fraction-to-boundary step length is not in (0, 1]");
      return step_length;
   }

   double PrimalDualInteriorPointMethod::dual_fraction_to_boundary(const OptimizationProblem& problem, const Multipliers& current_multipliers,
         Multipliers& direction_multipliers, double tau) {
      double step_length = 1.;
      for (const size_t variable_index: problem.get_lower_bounded_variables()) {
         if (direction_multipliers.lower_bounds[variable_index] < 0.) {
            double distance = -tau * current_multipliers.lower_bounds[variable_index] / direction_multipliers.lower_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      for (const size_t variable_index: problem.get_upper_bounded_variables()) {
         if (0. < direction_multipliers.upper_bounds[variable_index]) {
            double distance = -tau * current_multipliers.upper_bounds[variable_index] / direction_multipliers.upper_bounds[variable_index];
            if (0. < distance) {
               step_length = std::min(step_length, distance);
            }
         }
      }
      assert(0. < step_length && step_length <= 1. && "The dual fraction-to-boundary step length is not in (0, 1]");
      return step_length;
   }

   // generate the right-hand side
   void PrimalDualInteriorPointMethod::assemble_augmented_rhs(const Multipliers& current_multipliers, size_t number_variables, size_t number_constraints) {
      this->augmented_system.rhs.fill(0.);

      // objective gradient
      for (const auto [variable_index, derivative]: this->objective_gradient) {
         this->augmented_system.rhs[variable_index] -= derivative;
      }

      // constraint: evaluations and gradients
      for (size_t constraint_index: Range(number_constraints)) {
         // Lagrangian
         if (current_multipliers.constraints[constraint_index] != 0.) {
            for (const auto [variable_index, derivative]: this->constraint_jacobian[constraint_index]) {
               this->augmented_system.rhs[variable_index] += current_multipliers.constraints[constraint_index] * derivative;
            }
         }
         // constraints
         this->augmented_system.rhs[number_variables + constraint_index] = -this->constraints[constraint_index];
      }
      DEBUG2 << "RHS: "; print_vector(DEBUG2, view(this->augmented_system.rhs, 0, number_variables + number_constraints)); DEBUG << '\n';
   }

   void PrimalDualInteriorPointMethod::assemble_primal_dual_direction(const OptimizationProblem& barrier_problem, const Vector<double>& current_primals,
         const Multipliers& current_multipliers, Vector<double>& direction_primals, Multipliers& direction_multipliers) {
      // form the primal-dual direction
      direction_primals = view(this->augmented_system.solution, 0, barrier_problem.number_variables);
      // retrieve the duals with correct signs (note the minus sign)
      direction_multipliers.constraints = view(-this->augmented_system.solution, barrier_problem.number_variables,
         barrier_problem.number_variables + barrier_problem.number_constraints);
      this->compute_bound_dual_direction(barrier_problem, current_primals, current_multipliers, direction_primals, direction_multipliers);

      // "fraction-to-boundary" rule for primal variables and constraints multipliers
      const double tau = std::max(this->parameters.tau_min, 1. - this->barrier_parameter());
      const double primal_step_length = PrimalDualInteriorPointMethod::primal_fraction_to_boundary(barrier_problem, current_primals, direction_primals, tau);
      const double bound_dual_step_length = PrimalDualInteriorPointMethod::dual_fraction_to_boundary(barrier_problem, current_multipliers, direction_multipliers, tau);
      DEBUG << "Fraction-to-boundary rules:\n";
      DEBUG << "primal step length = " << primal_step_length << '\n';
      DEBUG << "bound dual step length = " << bound_dual_step_length << "\n\n";
      // scale the primal-dual variables
      direction_primals.scale(primal_step_length);
      direction_multipliers.constraints.scale(primal_step_length);
      direction_multipliers.lower_bounds.scale(bound_dual_step_length);
      direction_multipliers.upper_bounds.scale(bound_dual_step_length);

      // determine if the direction is a "small direction" (Section 3.9 of the Ipopt paper) TODO
      const bool is_small_step = PrimalDualInteriorPointMethod::is_small_step(barrier_problem, current_primals, direction_primals);
      if (is_small_step) {
         DEBUG << "This is a small step\n";
      }
   }

   void PrimalDualInteriorPointMethod::compute_bound_dual_direction(const OptimizationProblem& problem, const Vector<double>& current_primals,
         const Multipliers& current_multipliers, const Vector<double>& primal_direction, Multipliers& direction_multipliers) {
      direction_multipliers.lower_bounds.fill(0.);
      direction_multipliers.upper_bounds.fill(0.);
      for (const size_t variable_index: problem.get_lower_bounded_variables()) {
         const double distance_to_bound = current_primals[variable_index] - problem.variable_lower_bound(variable_index);
         direction_multipliers.lower_bounds[variable_index] = (this->barrier_parameter() - primal_direction[variable_index] *
            current_multipliers.lower_bounds[variable_index]) / distance_to_bound - current_multipliers.lower_bounds[variable_index];
         assert(is_finite(direction_multipliers.lower_bounds[variable_index]) && "The lower bound dual is infinite");
      }
      for (const size_t variable_index: problem.get_upper_bounded_variables()) {
         const double distance_to_bound = current_primals[variable_index] - problem.variable_upper_bound(variable_index);
         direction_multipliers.upper_bounds[variable_index] = (this->barrier_parameter() - primal_direction[variable_index] *
            current_multipliers.upper_bounds[variable_index]) / distance_to_bound - current_multipliers.upper_bounds[variable_index];
         assert(is_finite(direction_multipliers.upper_bounds[variable_index]) && "The upper bound dual is infinite");
      }
   }

   void PrimalDualInteriorPointMethod::compute_least_square_multipliers(const OptimizationProblem& problem, Iterate& iterate,
         Vector<double>& constraint_multipliers) {
      this->augmented_system.matrix.set_dimension(problem.number_variables + problem.number_constraints);
      this->augmented_system.matrix.reset();
      Preprocessing::compute_least_square_multipliers(problem.model, this->augmented_system.matrix, this->augmented_system.rhs, *this->linear_solver,
         iterate, constraint_multipliers, this->least_square_multiplier_max_norm);
   }

   void PrimalDualInteriorPointMethod::postprocess_iterate(const OptimizationProblem& problem, Vector<double>& primals, Multipliers& multipliers) {
      // rescale the bound multipliers (Eq. 16 in Ipopt paper)
      for (const size_t variable_index: problem.get_lower_bounded_variables()) {
         const double coefficient = this->barrier_parameter() / (primals[variable_index] - problem.variable_lower_bound(variable_index));
         if (is_finite(coefficient)) {
            const double lb = coefficient / this->parameters.k_sigma;
            const double ub = coefficient * this->parameters.k_sigma;
            if (lb <= ub) {
               const double current_value = multipliers.lower_bounds[variable_index];
               multipliers.lower_bounds[variable_index] = std::max(std::min(multipliers.lower_bounds[variable_index], ub), lb);
               if (multipliers.lower_bounds[variable_index] != current_value) {
                  DEBUG << "Multiplier for lower bound " << variable_index << " rescaled from " << current_value << " to " << multipliers.lower_bounds[variable_index] << '\n';
               }
            }
            else {
               WARNING << "Barrier subproblem: the bounds are in the wrong order in the lower bound multiplier reset\n";
            }
         }

      }
      for (const size_t variable_index: problem.get_upper_bounded_variables()) {
         const double coefficient = this->barrier_parameter() / (primals[variable_index] - problem.variable_upper_bound(variable_index));
         if (is_finite(coefficient)) {
            const double lb = coefficient * this->parameters.k_sigma;
            const double ub = coefficient / this->parameters.k_sigma;
            if (lb <= ub) {
               const double current_value = multipliers.upper_bounds[variable_index];
               multipliers.upper_bounds[variable_index] = std::max(std::min(multipliers.upper_bounds[variable_index], ub), lb);
               if (multipliers.upper_bounds[variable_index] != current_value) {
                  DEBUG << "Multiplier for upper bound " << variable_index << " rescaled from " << current_value << " to " << multipliers.upper_bounds[variable_index] << '\n';
               }
            }
            else {
               WARNING << "Barrier subproblem: the bounds are in the wrong order in the upper bound multiplier reset\n";
            }
         }
      }
   }

   void PrimalDualInteriorPointMethod::set_initial_point(const Vector<double>& /*point*/) {
      // do nothing
   }

   std::string PrimalDualInteriorPointMethod::get_name() const {
      return "primal-dual interior-point method";
   }
} // namespace
