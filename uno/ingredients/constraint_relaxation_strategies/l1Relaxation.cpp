// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/subproblems/Direction.hpp"
#include "ingredients/subproblems/Subproblem.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/VectorView.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

/*
 * Infeasibility detection and SQP methods for nonlinear optimization
 * Richard H. Byrd, Frank E. Curtis and Jorge Nocedal
 * http://epubs.siam.org/doi/pdf/10.1137/080738222
 */

namespace uno {
   l1Relaxation::l1Relaxation(const Model& model, const Options& options) :
         // call delegating constructor
         l1Relaxation(model,
               // create the l1 feasibility problem (objective multiplier = 0)
               l1RelaxedProblem(model, 0., options.get_double("l1_constraint_violation_coefficient"), 0., nullptr),
               // create the l1 relaxed problem
               l1RelaxedProblem(model, options.get_double("l1_relaxation_initial_parameter"), options.get_double("l1_constraint_violation_coefficient"),
                  0., nullptr),
               options) {
   }

   // private delegating constructor
   l1Relaxation::l1Relaxation(const Model& model, l1RelaxedProblem&& feasibility_problem, l1RelaxedProblem&& l1_relaxed_problem, const Options& options) :
         ConstraintRelaxationStrategy(model, l1_relaxed_problem.number_variables, l1_relaxed_problem.number_constraints,
               l1_relaxed_problem.number_objective_gradient_nonzeros(), l1_relaxed_problem.number_jacobian_nonzeros(),
               l1_relaxed_problem.number_hessian_nonzeros(), options),
         feasibility_problem(std::forward<l1RelaxedProblem>(feasibility_problem)),
         l1_relaxed_problem(std::forward<l1RelaxedProblem>(l1_relaxed_problem)),
         penalty_parameter(options.get_double("l1_relaxation_initial_parameter")),
         tolerance(options.get_double("tolerance")),
         parameters({
               options.get_bool("l1_relaxation_fixed_parameter"),
               options.get_double("l1_relaxation_decrease_factor"),
               options.get_double("l1_relaxation_epsilon1"),
               options.get_double("l1_relaxation_epsilon2"),
               options.get_double("l1_relaxation_residual_small_threshold")
         }),
         small_duals_threshold(options.get_double("l1_small_duals_threshold")),
         trial_multipliers(this->l1_relaxed_problem.number_variables, model.number_constraints) {
   }

   void l1Relaxation::initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) {
      // statistics
      this->subproblem->initialize_statistics(statistics, options);
      statistics.add_column("penalty param.", Statistics::double_width, options.get_int("statistics_penalty_parameter_column_order"));
      statistics.set("penalty param.", this->penalty_parameter);

      // initial iterate
      initial_iterate.feasibility_residuals.lagrangian_gradient.resize(this->feasibility_problem.number_variables);
      initial_iterate.feasibility_multipliers.lower_bounds.resize(this->feasibility_problem.number_variables);
      initial_iterate.feasibility_multipliers.upper_bounds.resize(this->feasibility_problem.number_variables);
      this->subproblem->set_elastic_variable_values(this->l1_relaxed_problem, initial_iterate);
      this->subproblem->generate_initial_iterate(this->l1_relaxed_problem, initial_iterate);
      this->evaluate_progress_measures(initial_iterate);
      this->compute_primal_dual_residuals(initial_iterate);
      this->set_statistics(statistics, initial_iterate);
      this->globalization_strategy->initialize(statistics, initial_iterate, options);
   }

   void l1Relaxation::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         WarmstartInformation& warmstart_information) {
      statistics.set("penalty param.", this->penalty_parameter);
      direction.reset();
      this->solve_sequence_of_relaxed_subproblems(statistics, current_iterate, direction, warmstart_information);
   }

   bool l1Relaxation::solving_feasibility_problem() const {
      return (this->penalty_parameter == 0.);
   }

   void l1Relaxation::switch_to_feasibility_problem(Statistics& /*statistics*/, Iterate& /*current_iterate*/) {
      throw std::runtime_error("l1Relaxation::switch_to_feasibility_problem is not implemented");
   }

   // use Byrd's steering rules to update the penalty parameter and compute a descent direction
   void l1Relaxation::solve_sequence_of_relaxed_subproblems(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         WarmstartInformation& warmstart_information) {
      // stage a: compute a direction for the current penalty parameter
      this->solve_l1_relaxed_problem(statistics, current_iterate, direction, this->penalty_parameter, warmstart_information);
      // from now on, only the penalty parameter, therefore the objective, changes
      warmstart_information.only_objective_changed();

      // penalty update: if penalty parameter is already 0 or fixed by the user, no need to decrease it
      if (0. < this->penalty_parameter && not this->parameters.fixed_parameter) {
         double linearized_residual = this->model.constraint_violation(current_iterate.evaluations.constraints +
               current_iterate.evaluations.constraint_jacobian * direction.primals, Norm::L1);
         DEBUG << "Linearized infeasibility mk(dk): " << linearized_residual << "\n\n";

         // if the current direction is already feasible, terminate
         if (this->tolerance < linearized_residual) {
            const double current_penalty_parameter = this->penalty_parameter;

            // stage c: compute the lowest possible constraint violation (penalty parameter = 0)
            DEBUG << "Compute ideal solution by solving the feasibility problem:\n";
            this->subproblem->initialize_feasibility_problem(this->feasibility_problem, current_iterate);
            Direction feasibility_direction(direction.number_variables, direction.number_constraints);
            this->solve_subproblem(statistics, this->feasibility_problem, current_iterate, current_iterate.feasibility_multipliers, feasibility_direction,
                  warmstart_information);
            std::swap(direction.multipliers, direction.feasibility_multipliers);
            const double residual_lowest_violation = this->model.constraint_violation(current_iterate.evaluations.constraints +
                  current_iterate.evaluations.constraint_jacobian * feasibility_direction.primals, Norm::L1);
            DEBUG << "Lowest linearized infeasibility mk(dk): " << residual_lowest_violation << '\n';
            this->subproblem->exit_feasibility_problem(this->feasibility_problem, current_iterate);

            // stage f: update the penalty parameter based on the current dual error
            this->decrease_parameter_aggressively(current_iterate, feasibility_direction);
            if (this->penalty_parameter < current_penalty_parameter) {
               this->solve_l1_relaxed_problem(statistics, current_iterate, direction, this->penalty_parameter, warmstart_information);
               linearized_residual = this->model.constraint_violation(current_iterate.evaluations.constraints +
                     current_iterate.evaluations.constraint_jacobian * direction.primals, Norm::L1);
            }

            // stage d: further decrease penalty parameter to reach a fraction of the ideal decrease
            this->enforce_linearized_residual_sufficient_decrease(statistics, current_iterate, direction, linearized_residual,
                  residual_lowest_violation, warmstart_information);
            // stage e: further decrease penalty parameter to guarantee a descent direction for the l1 merit function
            this->enforce_descent_direction_for_l1_merit(statistics, current_iterate, direction, feasibility_direction, warmstart_information);

            // save the dual feasibility direction
            direction.feasibility_multipliers = feasibility_direction.multipliers;
         }
      }
   }

   void l1Relaxation::solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, const WarmstartInformation& warmstart_information) {
      DEBUG << "Solving the subproblem with penalty parameter " << problem.get_objective_multiplier() << "\n\n";

      // solve the subproblem
      direction.set_dimensions(problem.number_variables, problem.number_constraints);
      this->subproblem->solve(statistics, problem, current_iterate, current_multipliers, direction, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, this->model.number_variables));
      DEBUG3 << direction << '\n';
      assert(direction.status == SubproblemStatus::OPTIMAL && "The subproblem was not solved to optimality");
   }

   void l1Relaxation::solve_l1_relaxed_problem(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double current_penalty_parameter, const WarmstartInformation& warmstart_information) {
      this->l1_relaxed_problem.set_objective_multiplier(current_penalty_parameter);
      this->solve_subproblem(statistics, this->l1_relaxed_problem, current_iterate, current_iterate.multipliers, direction, warmstart_information);
   }

   void l1Relaxation::decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction) {
      this->trial_multipliers.constraints = current_iterate.feasibility_multipliers.constraints + direction.feasibility_multipliers.constraints;
      this->trial_multipliers.lower_bounds = current_iterate.feasibility_multipliers.lower_bounds + direction.feasibility_multipliers.lower_bounds;
      this->trial_multipliers.upper_bounds = current_iterate.feasibility_multipliers.upper_bounds + direction.feasibility_multipliers.upper_bounds;

      // there must be at least a nonzero dual to avoid trivial stationary points
      if (this->trial_multipliers.not_all_zero(this->model.number_variables, this->small_duals_threshold)) {
         // compute the ideal error (with a zero penalty parameter)
         const double infeasible_dual_error = l1Relaxation::compute_infeasible_dual_error(current_iterate);
         DEBUG << "Ideal dual error: " << infeasible_dual_error << '\n';
         const double scaled_error = infeasible_dual_error / std::max(1., current_iterate.primal_feasibility);
         this->penalty_parameter = std::min(this->penalty_parameter, scaled_error * scaled_error);
         DEBUG << "Further aggressively decrease the penalty parameter to " << this->penalty_parameter << '\n';
      }
      else {
         DEBUG << RED << "l1Relaxation: all multipliers are almost 0. The penalty parameter won't be decreased" << RESET << '\n';
      }
   }

   // measure that combines KKT error and complementarity error
   double l1Relaxation::compute_infeasible_dual_error(Iterate& current_iterate) {
      // stationarity error
      this->feasibility_problem.evaluate_lagrangian_gradient(current_iterate.feasibility_residuals.lagrangian_gradient, current_iterate, this->trial_multipliers);
      double error = norm_1(current_iterate.residuals.lagrangian_gradient.constraints_contribution);

      // complementarity error
      const double shift_value = 0.;
      error += this->feasibility_problem.complementarity_error(current_iterate.primals, current_iterate.evaluations.constraints,
            this->trial_multipliers, shift_value, Norm::L1);
      return error;
   }

   void l1Relaxation::enforce_linearized_residual_sufficient_decrease(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double linearized_residual, double residual_lowest_violation, WarmstartInformation& warmstart_information) {
      while (0. < this->penalty_parameter && not this->linearized_residual_sufficient_decrease(current_iterate, linearized_residual,
            residual_lowest_violation)) {
         // decrease the penalty parameter and re-solve the problem
         this->penalty_parameter /= this->parameters.decrease_factor;
         DEBUG << "Further decrease the penalty parameter to " << this->penalty_parameter << '\n';
         this->solve_l1_relaxed_problem(statistics, current_iterate, direction, this->penalty_parameter, warmstart_information);

         // recompute the linearized residual
         linearized_residual = this->model.constraint_violation(current_iterate.evaluations.constraints +
               current_iterate.evaluations.constraint_jacobian * direction.primals, Norm::L1);
         DEBUG << "Linearized infeasibility mk(dk): " << linearized_residual << "\n\n";
      }
      DEBUG << "Condition enforce_linearized_residual_sufficient_decrease is true\n";
   }

   bool l1Relaxation::linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual,
         double residual_lowest_violation) const {
      if (residual_lowest_violation <= this->parameters.residual_small_threshold) {
         return (linearized_residual <= this->parameters.residual_small_threshold);
      }
      const double linearized_residual_reduction = current_iterate.progress.infeasibility - linearized_residual;
      const double lowest_linearized_residual_reduction = current_iterate.progress.infeasibility - residual_lowest_violation;
      if (lowest_linearized_residual_reduction < 0.) {
         WARNING << "The lowest_linearized_residual_reduction quantity is negative. Negative curvature in your problem\n";
      }
      return (linearized_residual_reduction >= this->parameters.epsilon1 * lowest_linearized_residual_reduction);
   }

   void l1Relaxation::enforce_descent_direction_for_l1_merit(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         const Direction& feasibility_direction, WarmstartInformation& warmstart_information) {
      while (0. < this->penalty_parameter && not this->is_descent_direction_for_l1_merit_function(current_iterate, direction, feasibility_direction)) {
         // decrease the penalty parameter and re-solve the problem
         this->penalty_parameter /= this->parameters.decrease_factor;
         DEBUG << "Further decrease the penalty parameter to " << this->penalty_parameter << '\n';
         this->solve_l1_relaxed_problem(statistics, current_iterate, direction, this->penalty_parameter, warmstart_information);
      }
      DEBUG << "Condition enforce_descent_direction_for_l1_merit is true\n\n";
   }

   bool l1Relaxation::is_descent_direction_for_l1_merit_function(const Iterate& current_iterate, const Direction& direction,
         const Direction& feasibility_direction) const {
      const double predicted_l1_merit_reduction = current_iterate.primal_feasibility - direction.subproblem_objective;
      const double lowest_decrease_objective = current_iterate.primal_feasibility - feasibility_direction.subproblem_objective;
      return (predicted_l1_merit_reduction >= this->parameters.epsilon2 * lowest_decrease_objective);
   }

   bool l1Relaxation::is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) {
      this->subproblem->postprocess_iterate(this->l1_relaxed_problem, trial_iterate);
      this->compute_progress_measures(current_iterate, trial_iterate);
      trial_iterate.objective_multiplier = this->l1_relaxed_problem.get_objective_multiplier();

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_iterate.evaluate_objective(this->model);
         accept_iterate = true;
         statistics.set("status", "0 primal step");
      }
      else {
         // invoke the globalization strategy for acceptance
         const ProgressMeasures predicted_reduction = this->compute_predicted_reduction_models(current_iterate, direction, step_length);
         accept_iterate = this->globalization_strategy->is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
               predicted_reduction, this->penalty_parameter);
      }
      if (accept_iterate) {
         this->check_exact_relaxation(trial_iterate);
         // this->set_dual_residuals_statistics(statistics, trial_iterate);
      }
      this->set_progress_statistics(statistics, trial_iterate);
      return accept_iterate;
   }

   void l1Relaxation::compute_primal_dual_residuals(Iterate& iterate) {
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(this->l1_relaxed_problem, this->feasibility_problem, iterate);
   }

   void l1Relaxation::evaluate_progress_measures(Iterate& iterate) const {
      this->set_infeasibility_measure(iterate);
      this->set_objective_measure(iterate);
      this->subproblem->set_auxiliary_measure(this->model, iterate);
   }

   ProgressMeasures l1Relaxation::compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length) {
      return {
         this->compute_predicted_infeasibility_reduction_model(current_iterate, direction.primals, step_length),
         this->first_order_predicted_reduction ? this->compute_predicted_objective_reduction_model(current_iterate, direction.primals, step_length) :
            this->compute_predicted_objective_reduction_model(current_iterate, direction.primals, step_length, this->subproblem->get_lagrangian_hessian()),
         this->subproblem->compute_predicted_auxiliary_reduction_model(this->model, current_iterate, direction.primals, step_length)
      };
   }

   size_t l1Relaxation::maximum_number_variables() const {
      return this->l1_relaxed_problem.number_variables;
   }

   size_t l1Relaxation::maximum_number_constraints() const {
      return this->l1_relaxed_problem.number_constraints;
   }

   // for information purposes, check that l1 is an exact relaxation
   void l1Relaxation::check_exact_relaxation(Iterate& iterate) const {
      const double norm_inf_multipliers = norm_inf(iterate.multipliers.constraints);
      if (0. < norm_inf_multipliers && this->penalty_parameter <= 1./norm_inf_multipliers) {
         DEBUG << "The value of the penalty parameter is consistent with an exact relaxation\n\n";
      }
   }

   void l1Relaxation::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const {
      statistics.set("stationarity", iterate.residuals.stationarity);
      statistics.set("complementarity", iterate.residuals.complementarity);
   }
} // namespace
