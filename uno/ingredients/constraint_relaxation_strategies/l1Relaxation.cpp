// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategyFactory.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationSpace.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

/*
 * A Sequential Quadratic Optimization Algorithm with Rapid Infeasibility Detection
 * James V. Burke, Frank E. Curtis, and Hao Wang
 * https://epubs.siam.org/doi/abs/10.1137/120880045
 */

namespace uno {
   l1Relaxation::l1Relaxation(const Options& options) :
         ConstraintRelaxationStrategy(options),
         inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         inertia_correction_strategy(InertiaCorrectionStrategyFactory::create(options)),
         penalty_parameter(0.1), // TODO add option
         tolerance(options.get_double("primal_tolerance")),
         parameters({
            1e-2, /* beta */
            0.1, /* theta */
            10., /* kappa_rho */
            10., /* kappa_lambda */
            1e-2, /* epsilon */
            1. - 1e-18, /* omega */
            0.5 /* delta */
         }) {
   }

   void l1Relaxation::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate,
         Direction& direction, double trust_region_radius, const Options& options) {
      this->relaxed_problem = std::make_unique<l1RelaxedProblem>(model, this->penalty_parameter, 1.);
      assert(this->relaxed_problem != nullptr);
      this->feasibility_problem = std::make_unique<const l1RelaxedProblem>(model, 0., 1.);
      assert(this->feasibility_problem != nullptr);

      // Hessian model
      this->hessian_model = HessianModelFactory::create(model, options);

      // memory allocation
      this->inequality_handling_method->initialize(*this->relaxed_problem, initial_iterate, *this->hessian_model,
         *this->inertia_correction_strategy, trust_region_radius);
      direction = Direction(this->relaxed_problem->number_variables, this->relaxed_problem->number_constraints);

      // statistics
      this->inertia_correction_strategy->initialize_statistics(statistics, options);
      this->inequality_handling_method->initialize_statistics(statistics, options);
      statistics.add_column("penalty param.", Statistics::double_width, options.get_int("statistics_penalty_parameter_column_order"));
      statistics.set("penalty param.", this->penalty_parameter);

      // initial iterate
      initial_iterate.set_number_variables(this->feasibility_problem->number_variables);
      initial_iterate.residuals.lagrangian_gradient.resize(this->feasibility_problem->number_variables);
      initial_iterate.multipliers.lower_bounds.resize(this->feasibility_problem->number_variables);
      initial_iterate.multipliers.upper_bounds.resize(this->feasibility_problem->number_variables);
      this->inequality_handling_method->generate_initial_iterate(initial_iterate);
      this->inequality_handling_method->initialize_feasibility_problem(initial_iterate);
      this->inequality_handling_method->set_elastic_variable_values(*this->feasibility_problem, initial_iterate);
      // initialize the feasibility multipliers
      this->feasibility_multipliers.constraints.resize(this->feasibility_problem->number_constraints);
      this->feasibility_multipliers.lower_bounds.resize(this->feasibility_problem->number_variables);
      this->feasibility_multipliers.upper_bounds.resize(this->feasibility_problem->number_variables);

      initial_iterate.evaluate_objective_gradient(model);
      initial_iterate.evaluate_constraints(model);
      this->inequality_handling_method->evaluate_constraint_jacobian(initial_iterate);
      this->relaxed_problem->evaluate_lagrangian_gradient(initial_iterate.residuals.lagrangian_gradient, *this->inequality_handling_method,
         initial_iterate);
      this->compute_primal_dual_residuals(*this->relaxed_problem, initial_iterate);
   }

   void l1Relaxation::compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& /*globalization_strategy*/,
         Iterate& current_iterate, Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information) {
      statistics.set("penalty param.", this->penalty_parameter);
      direction.reset();

      // solve feasibility problem
      DEBUG << "Solving the l1 feasibility problem\n";
      Direction feasibility_direction(this->feasibility_problem->number_variables, this->feasibility_problem->number_constraints);
      std::swap(current_iterate.multipliers, this->feasibility_multipliers);
      this->feasibility_inequality_handling_method->solve(statistics, current_iterate, feasibility_direction, *this->hessian_model,
         *this->inertia_correction_strategy, trust_region_radius, warmstart_information);
      assert(direction.status == SubproblemStatus::OPTIMAL && "The feasibility subproblem was not solved to optimality");
      std::swap(current_iterate.multipliers, this->feasibility_multipliers);
      // assemble multipliers for feasibility problem
      Multipliers trial_feasibility_multipliers(current_iterate.number_variables, current_iterate.number_constraints);
      trial_feasibility_multipliers.constraints = this->feasibility_multipliers.constraints + feasibility_direction.multipliers.constraints;
      trial_feasibility_multipliers.lower_bounds = this->feasibility_multipliers.lower_bounds + feasibility_direction.multipliers.lower_bounds;
      trial_feasibility_multipliers.upper_bounds = this->feasibility_multipliers.upper_bounds + feasibility_direction.multipliers.upper_bounds;
      const double feasibility_dual_error = this->Rinf(current_iterate, trial_feasibility_multipliers);
      DEBUG << "Infeasibility dual error = " << feasibility_dual_error << '\n';

      // update penalty parameter and duals
      // test equation (3.12)
      if (this->tolerance < this->v(current_iterate) &&
            this->delta_l(feasibility_direction, current_iterate) <= this->parameters.theta * this->v(current_iterate)) {
         // update penalty parameter according to (3.13)
         this->penalty_parameter = std::min(this->penalty_parameter, this->parameters.kappa_rho * feasibility_dual_error * feasibility_dual_error);
         DEBUG << "Penalty parameter updated to " << this->penalty_parameter << '\n';

         // update current multipliers according to (3.14)
         const double multipliers_distance = norm_2(
               current_iterate.multipliers.constraints + (-1) * this->feasibility_multipliers.constraints);
         const double alpha = std::min(1., this->parameters.kappa_lambda * feasibility_dual_error * feasibility_dual_error / multipliers_distance);
         current_iterate.multipliers.constraints = alpha * current_iterate.multipliers.constraints +
            (1. - alpha) * this->feasibility_multipliers.constraints;
         DEBUG2 << "Updated multipliers: " << current_iterate.multipliers.constraints << '\n';
      }

      // solve l1 relaxed problem
      DEBUG << "\nSolving the regular l1 relaxed problem\n";
      this->relaxed_problem->set_objective_multiplier(this->penalty_parameter);
      //this->subproblem->set_initial_point(feasibility_direction.primals);
      warmstart_information.only_objective_changed();
      Direction optimality_direction(direction.number_variables, direction.number_constraints);
      this->inequality_handling_method->solve(statistics, current_iterate, optimality_direction, *this->hessian_model,
         *this->inertia_correction_strategy, trust_region_radius, warmstart_information);
      assert(direction.status == SubproblemStatus::OPTIMAL && "The l1 relaxed subproblem was not solved to optimality");
      // assemble multipliers for l1 relaxed problem
      Multipliers trial_multipliers(current_iterate.number_variables, current_iterate.number_constraints);
      trial_multipliers.constraints = current_iterate.multipliers.constraints + optimality_direction.multipliers.constraints;
      trial_multipliers.lower_bounds = current_iterate.multipliers.lower_bounds + optimality_direction.multipliers.lower_bounds;
      trial_multipliers.upper_bounds = current_iterate.multipliers.upper_bounds + optimality_direction.multipliers.upper_bounds;
      const double optimality_dual_error = this->Ropt(current_iterate, this->penalty_parameter, trial_multipliers);
      DEBUG << "Optimality dual error = " << optimality_dual_error << '\n';

      // interpolate between two directions
      DEBUG << "\nInterpolating between the two directions:\n";
      DEBUG2 << "Feasibility direction: ";
      print_vector(DEBUG2, view(feasibility_direction.primals, 0, this->relaxed_problem->model.number_variables));
      DEBUG2 << "Optimality direction:  ";
      print_vector(DEBUG2, view(optimality_direction.primals, 0, this->relaxed_problem->model.number_variables));
      const double w = this->compute_w(feasibility_direction, optimality_direction, current_iterate);
      const auto interpolated_direction = w * feasibility_direction.primals + (1 - w) * optimality_direction.primals;
      for (size_t variable_index: Range(direction.number_variables)) {
         direction.primals[variable_index] = interpolated_direction[variable_index];
      }

      // update the penalty parameter (3.18)
      const double multipliers_inf_norm = norm_inf(trial_multipliers.constraints, trial_multipliers.lower_bounds, trial_multipliers.upper_bounds);
      if (this->penalty_parameter * multipliers_inf_norm > 1.) {
         this->penalty_parameter = std::min(this->parameters.delta * this->penalty_parameter, (1. - this->parameters.epsilon) / multipliers_inf_norm);
         DEBUG << "Penalty parameter updated to " << this->penalty_parameter << '\n';
      }
      // update the penalty parameter (3.19)
      const double delta_m = this->delta_m(direction, current_iterate, this->penalty_parameter);
      const double delta_l = this->delta_l(direction, current_iterate);
      if (delta_m >= this->parameters.epsilon * delta_l && w >= this->parameters.omega) {
         this->penalty_parameter *= this->parameters.delta;
      }
      else if (delta_m < this->parameters.epsilon * delta_l) {
         const double zeta = this->compute_zeta(direction, current_iterate);
         this->penalty_parameter = std::min(this->parameters.delta * this->penalty_parameter, zeta);
      }
      DEBUG << "Penalty parameter updated to " << this->penalty_parameter << '\n';

      // construct the final direction
      direction.multipliers = optimality_direction.multipliers;
      direction.norm = norm_inf(view(direction.primals, 0, this->relaxed_problem->model.number_variables));
      DEBUG3 << direction << '\n';
      warmstart_information.no_changes();
   }

   bool l1Relaxation::solving_feasibility_problem() const {
      return (this->penalty_parameter == 0.);
   }

   void l1Relaxation::switch_to_feasibility_problem(Statistics& /*statistics*/, GlobalizationStrategy& /*globalization_strategy*/,
         Iterate& /*current_iterate*/, double /*trust_region_radius*/,
         WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("l1Relaxation::switch_to_feasibility_problem is not implemented\n");
   }

   bool l1Relaxation::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         double trust_region_radius, const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      const bool accept_iterate = this->inequality_handling_method->is_iterate_acceptable(statistics, globalization_strategy,
         *this->hessian_model, *this->inertia_correction_strategy, trust_region_radius, current_iterate, trial_iterate,
         direction, step_length, user_callbacks);
      if (accept_iterate) {
         this->check_exact_relaxation(trial_iterate);
      }
      trial_iterate.status = this->check_termination(model, trial_iterate);
      warmstart_information.no_changes();
      return accept_iterate;
   }

   SolutionStatus l1Relaxation::check_termination(const Model& model, Iterate& iterate) {
      iterate.evaluate_objective_gradient(model);
      iterate.evaluate_constraints(model);

      this->relaxed_problem->evaluate_lagrangian_gradient(iterate.residuals.lagrangian_gradient, *this->inequality_handling_method,
         iterate);
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(*this->relaxed_problem, iterate);
      return ConstraintRelaxationStrategy::check_termination(*this->relaxed_problem, iterate);
   }

   std::string l1Relaxation::get_name() const {
      return this->inequality_handling_method->get_name() + " with " + this->hessian_model->name + " Hessian and " +
         this->inertia_correction_strategy->get_name() + " regularization";
   }

   size_t l1Relaxation::get_number_subproblems_solved() const {
      return this->inequality_handling_method->number_subproblems_solved;
   }

   // protected member functions

   // v: l1 norm: ||c(x_k)||_1
   double l1Relaxation::v(const Iterate& current_iterate) const {
      return current_iterate.progress.infeasibility;
   }

   // l: constraint violation of linearized constraint: ||c(x_k) + \nabla c(x_k)^T d||_1
   double l1Relaxation::l(const Direction& direction, const Iterate& current_iterate) const {
      // TODO preallocate
      Vector<double> constraint_directional_derivative(this->relaxed_problem->model.number_constraints);
      const auto& evaluation_space = this->inequality_handling_method->get_evaluation_space();
      evaluation_space.compute_constraint_jacobian_vector_product(direction.primals, constraint_directional_derivative);
      const auto linearized_constraint = current_iterate.evaluations.constraints + constraint_directional_derivative;
      return this->relaxed_problem->model.constraint_violation(linearized_constraint, Norm::L1);
   }

   // delta_l: predicted (linear model) reduction of constraint violation
   // ||c(x_k)||_1 - ||c(x_k) + \alpha \nabla c(x_k)^T d||_1
   double l1Relaxation::delta_l(const Direction& direction, const Iterate& current_iterate) const {
      return this->v(current_iterate) - this->l(direction, current_iterate);
   }

   // Ropt: measure that combines stationarity error and complementarity error for l1 relaxed problem
   double l1Relaxation::Ropt(Iterate& current_iterate, double objective_multiplier, const Multipliers& multipliers) const {
      /*
      // stationarity error
      this->relaxed_problem.evaluate_lagrangian_gradient(current_iterate.residuals.lagrangian_gradient, current_iterate, multipliers);
      const auto scaled_lagrangian = objective_multiplier * current_iterate.residuals.lagrangian_gradient.objective_contribution +
         current_iterate.residuals.lagrangian_gradient.constraints_contribution;
      double error = norm_inf(scaled_lagrangian);

      // complementarity error
      error += this->relaxed_problem.complementarity_error(current_iterate.primals, current_iterate.evaluations.constraints, multipliers, 0., Norm::INF);
      return error;
      */
      return 0.;
   }

   // Rinf: measure that combines stationarity error and complementarity error for feasibility problem
   double l1Relaxation::Rinf(Iterate& current_iterate, const Multipliers& multipliers) const {
      /*
      // stationarity error
      this->feasibility_problem.evaluate_lagrangian_gradient(current_iterate.feasibility_residuals.lagrangian_gradient, current_iterate, multipliers);
      const auto scaled_lagrangian = current_iterate.feasibility_residuals.lagrangian_gradient.constraints_contribution;
      double error = norm_inf(scaled_lagrangian);

      // complementarity error
      error += this->feasibility_problem.complementarity_error(current_iterate.primals, current_iterate.evaluations.constraints, multipliers, 0., Norm::INF);
      return error;
      */
      return 0.;
   }

   // m: predicted reduction of the l1 merit function
   double l1Relaxation::delta_m(const Direction& direction, const Iterate& current_iterate, double objective_multiplier) const {
      return -objective_multiplier * dot(direction.primals, current_iterate.evaluations.objective_gradient) +
         this->delta_l(direction, current_iterate);
   }

   double l1Relaxation::compute_w(const Direction& feasibility_direction, const Direction& optimality_direction,
         const Iterate& current_iterate) {
      const double delta_l_dbar = this->delta_l(feasibility_direction, current_iterate);
      // TODO preallocate
      Vector<double> trial_direction(optimality_direction.primals.size());
      Vector<double> constraint_directional_derivative(this->relaxed_problem->model.number_constraints);
      const auto& evaluation_space = this->inequality_handling_method->get_evaluation_space();
      double weight = 0.;
      double lower_bound = 0.;
      double upper_bound = 1.;
      bool termination = false;
      while (!termination) {
         DEBUG2 << "Testing the interpolation weight " << weight << '\n';
         trial_direction = weight * feasibility_direction.primals + (1 - weight) * optimality_direction.primals;
         // update reduction in linearized feasibility model
         evaluation_space.compute_constraint_jacobian_vector_product(trial_direction, constraint_directional_derivative);
         const auto trial_linearized_constraints = current_iterate.evaluations.constraints + constraint_directional_derivative;
         double delta_l_d = current_iterate.progress.infeasibility - this->relaxed_problem->model.constraint_violation(trial_linearized_constraints, Norm::L1);
         DEBUG2 << "Trial predicted infeasibility reduction = " << delta_l_d << '\n';

         // sufficient decrease condition
         if (delta_l_d >= this->parameters.beta * delta_l_dbar) {
            // test if bisection has succeeded
            if (weight == 0. || upper_bound - lower_bound <= 1e-8) {
               termination = true;
               DEBUG << "Decrease condition satisfied, terminate with interpolation weight " << weight << '\n';
               DEBUG2 << "Interpolated direction:  ";
               print_vector(DEBUG2, trial_direction);
            }
            else {
               // keep the first half
               upper_bound = weight;
               weight = (lower_bound + upper_bound) / 2.;
               DEBUG2 << "Decrease condition satisfied, decreasing the interpolation weight\n";
            }
         }
         else if (weight == 1.) {
            termination = true;
            DEBUG << "Terminate with interpolation weight " << weight << '\n';
            DEBUG2 << "Interpolated direction:  ";
            print_vector(DEBUG2, trial_direction);
         }
         else {
            // sufficient decrease condition violated: keep the second half
            lower_bound = weight;
            weight = (lower_bound + upper_bound) / 2.;
            DEBUG2 << "Decrease condition violated, increasing the interpolation weight\n";
         }
      }
      return weight;
   }

   // zeta: upper bound on objective multiplier
   double l1Relaxation::compute_zeta(const Direction& direction, const Iterate& current_iterate) const {
      const auto& evaluation_space = this->inequality_handling_method->get_evaluation_space();
      const double numerator = (1. - this->parameters.epsilon) * this->delta_l(direction, current_iterate);
      const double denominator = dot(direction.primals, current_iterate.evaluations.objective_gradient);
         /*+ this->inequality_handling_method->compute_hessian_quadratic_product()
         evaluation_space.compute_hessian_quadratic_product(subproblem, direction.primals) / 2.;*/
      return numerator / denominator;
   }

   // for information purposes, check that l1 is an exact relaxation
   void l1Relaxation::check_exact_relaxation(const Iterate& iterate) const {
      const double norm_inf_multipliers = norm_inf(iterate.multipliers.constraints);
      if (0. < norm_inf_multipliers && this->penalty_parameter <= 1. / norm_inf_multipliers) {
         DEBUG << "The value of the penalty parameter is consistent with an exact relaxation\n";
      }
   }
} // namespace