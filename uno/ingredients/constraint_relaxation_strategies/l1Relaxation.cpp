// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/constraint_relaxation_strategies/l1RelaxedProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategyFactory.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/Expression.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

/*
 * Infeasibility detection and SQP methods for nonlinear optimization
 * Richard H. Byrd, Frank E. Curtis and Jorge Nocedal
 * http://epubs.siam.org/doi/pdf/10.1137/080738222
 */

namespace uno {
   l1Relaxation::l1Relaxation(const Options& options):
         ConstraintRelaxationStrategy(options),
         penalty_parameter(options.get_double("l1_relaxation_initial_parameter")),
         constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")),
         l1_relaxed_hessian_model(HessianModelFactory::create(options)),
         feasibility_hessian_model(HessianModelFactory::create(options)),
         l1_relaxed_regularization_strategy(RegularizationStrategyFactory::create(options)),
         feasibility_regularization_strategy(RegularizationStrategyFactory::create(options)),
         inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         feasibility_inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         tolerance(options.get_double("tolerance")),
         parameters({
               options.get_bool("l1_relaxation_fixed_parameter"),
               options.get_double("l1_relaxation_decrease_factor"),
               options.get_double("l1_relaxation_epsilon1"),
               options.get_double("l1_relaxation_epsilon2"),
               options.get_double("l1_relaxation_residual_small_threshold")
         }),
         small_duals_threshold(options.get_double("l1_small_duals_threshold")) {
   }

   void l1Relaxation::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, Direction& direction,
         const Options& options) {
      const l1RelaxedProblem l1_relaxed_problem{model, this->penalty_parameter, this->constraint_violation_coefficient};
      const l1RelaxedProblem feasibility_problem{model, 0., this->constraint_violation_coefficient};

      // memory allocation
      this->l1_relaxed_hessian_model->initialize(model);
      this->feasibility_hessian_model->initialize(model);
      this->inequality_handling_method->initialize(l1_relaxed_problem, *this->l1_relaxed_hessian_model,
         *this->l1_relaxed_regularization_strategy);
      this->feasibility_inequality_handling_method->initialize(feasibility_problem, *this->feasibility_hessian_model,
         *this->feasibility_regularization_strategy);
      this->trial_multipliers.constraints.resize(l1_relaxed_problem.number_constraints);
      this->trial_multipliers.lower_bounds.resize(l1_relaxed_problem.number_variables);
      this->trial_multipliers.upper_bounds.resize(l1_relaxed_problem.number_variables);
      direction = Direction(l1_relaxed_problem.number_variables, l1_relaxed_problem.number_constraints);

      // statistics
      this->l1_relaxed_regularization_strategy->initialize_statistics(statistics, options);
      this->feasibility_regularization_strategy->initialize_statistics(statistics, options);
      this->inequality_handling_method->initialize_statistics(statistics, options);
      this->feasibility_inequality_handling_method->initialize_statistics(statistics, options);
      statistics.add_column("penalty", Statistics::double_width - 5, options.get_int("statistics_penalty_parameter_column_order"));
      statistics.set("penalty", this->penalty_parameter);

      // initial iterate
      initial_iterate.feasibility_residuals.lagrangian_gradient.resize(l1_relaxed_problem.number_variables);
      initial_iterate.feasibility_multipliers.lower_bounds.resize(l1_relaxed_problem.number_variables);
      initial_iterate.feasibility_multipliers.upper_bounds.resize(l1_relaxed_problem.number_variables);
      this->inequality_handling_method->set_elastic_variable_values(l1_relaxed_problem, initial_iterate);
      this->inequality_handling_method->generate_initial_iterate(l1_relaxed_problem, initial_iterate);
      this->evaluate_progress_measures(*this->inequality_handling_method, l1_relaxed_problem, initial_iterate);
      this->compute_primal_dual_residuals(model, initial_iterate);
      this->set_statistics(statistics, model, initial_iterate);
   }

   void l1Relaxation::compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& /*globalization_strategy*/, const Model& model,
         Iterate& current_iterate, Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information) {
      statistics.set("penalty", this->penalty_parameter);
      direction.reset();
      this->solve_sequence_of_relaxed_subproblems(statistics, model, current_iterate, direction, trust_region_radius, warmstart_information);
   }
   
   bool l1Relaxation::solving_feasibility_problem() const {
      return (this->penalty_parameter == 0.);
   }

   void l1Relaxation::switch_to_feasibility_problem(Statistics& /*statistics*/, GlobalizationStrategy& /*globalization_strategy*/,
         const Model& /*model*/, Iterate& /*current_iterate*/, WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("l1Relaxation::switch_to_feasibility_problem is not implemented");
   }

   // use Byrd's steering rules to update the penalty parameter and compute a descent direction
   void l1Relaxation::solve_sequence_of_relaxed_subproblems(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information) {
      // stage a: compute a direction for the current penalty parameter
      this->solve_l1_relaxed_problem(statistics, model, current_iterate, direction, this->penalty_parameter, trust_region_radius,
         warmstart_information);

      // penalty update: if penalty parameter is already 0 or fixed by the user, no need to decrease it
      if (0. < this->penalty_parameter && !this->parameters.fixed_parameter) {
         double linearized_residual = model.constraint_violation(current_iterate.evaluations.constraints +
            current_iterate.evaluations.constraint_jacobian * direction.primals, Norm::L1);
         DEBUG << "Linearized infeasibility mk(dk): " << linearized_residual << "\n\n";

         // terminate if the current direction is already feasible, otherwise adjust the penalty parameter
         if (this->tolerance < linearized_residual) {
            const double current_penalty_parameter = this->penalty_parameter;

            // stage c: compute the lowest possible constraint violation (penalty parameter = 0)
            DEBUG << "Compute ideal solution by solving the feasibility problem:\n";
            const l1RelaxedProblem feasibility_problem{model, 0., this->constraint_violation_coefficient};
            this->feasibility_inequality_handling_method->initialize_feasibility_problem(feasibility_problem, current_iterate);
            Direction feasibility_direction(direction.number_variables, direction.number_constraints);
            this->solve_subproblem(statistics, *this->feasibility_inequality_handling_method, feasibility_problem, current_iterate,
               current_iterate.feasibility_multipliers, feasibility_direction, *this->feasibility_hessian_model,
               *this->feasibility_regularization_strategy, trust_region_radius, warmstart_information);
            std::swap(feasibility_direction.multipliers, feasibility_direction.feasibility_multipliers);
            const double residual_lowest_violation = model.constraint_violation(current_iterate.evaluations.constraints +
               current_iterate.evaluations.constraint_jacobian * feasibility_direction.primals, Norm::L1);
            DEBUG << "Lowest linearized infeasibility mk(dk): " << residual_lowest_violation << '\n';
            this->feasibility_inequality_handling_method->exit_feasibility_problem(feasibility_problem, current_iterate);

            // stage f: update the penalty parameter based on the current dual error
            this->decrease_parameter_aggressively(model, current_iterate, feasibility_direction);
            if (this->penalty_parameter < current_penalty_parameter) {
               this->solve_l1_relaxed_problem(statistics, model, current_iterate, direction, this->penalty_parameter,
                  trust_region_radius, warmstart_information);
               linearized_residual = model.constraint_violation(current_iterate.evaluations.constraints +
                  current_iterate.evaluations.constraint_jacobian * direction.primals, Norm::L1);
            }

            // stage d: further decrease penalty parameter to reach a fraction of the ideal decrease
            this->enforce_linearized_residual_sufficient_decrease(statistics, model, current_iterate, direction,
               linearized_residual, residual_lowest_violation, trust_region_radius, warmstart_information);
            // stage e: further decrease penalty parameter to guarantee a descent direction for the l1 merit function
            this->enforce_descent_direction_for_l1_merit(statistics, model, current_iterate, direction,
               feasibility_direction, trust_region_radius, warmstart_information);

            // save the dual feasibility direction
            direction.feasibility_multipliers = feasibility_direction.feasibility_multipliers;
         }
      }
   }

   void l1Relaxation::solve_subproblem(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers, Direction& direction,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      DEBUG << "Solving the subproblem with penalty parameter " << problem.get_objective_multiplier() << "\n\n";

      // solve the subproblem
      direction.set_dimensions(problem.number_variables, problem.number_constraints);
      inequality_handling_method.solve(statistics, problem, current_iterate, current_multipliers, direction, hessian_model,
         regularization_strategy, trust_region_radius, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
      assert(direction.status == SubproblemStatus::OPTIMAL && "The subproblem was not solved to optimality");
   }

   void l1Relaxation::solve_l1_relaxed_problem(Statistics& statistics, const Model& model, Iterate& current_iterate, Direction& direction,
         double current_penalty_parameter, double trust_region_radius, WarmstartInformation& warmstart_information) {
      const l1RelaxedProblem l1_relaxed_problem{model, current_penalty_parameter, this->constraint_violation_coefficient};
      this->solve_subproblem(statistics, *this->inequality_handling_method, l1_relaxed_problem, current_iterate,
         current_iterate.multipliers, direction, *this->l1_relaxed_hessian_model, *this->l1_relaxed_regularization_strategy,
         trust_region_radius, warmstart_information);
      if (direction.status == SubproblemStatus::UNBOUNDED_PROBLEM) {
         throw std::runtime_error("l1Relaxation::solve_l1_relaxed_problem: the subproblem is unbounded, this should not happen. "
            "If the subproblem has curvature, use regularization. If not, use a trust-region method.\n");
      }
   }

   void l1Relaxation::decrease_parameter_aggressively(const Model& model, Iterate& current_iterate, const Direction& direction) {
      this->trial_multipliers.constraints = current_iterate.feasibility_multipliers.constraints + direction.feasibility_multipliers.constraints;
      this->trial_multipliers.lower_bounds = current_iterate.feasibility_multipliers.lower_bounds + direction.feasibility_multipliers.lower_bounds;
      this->trial_multipliers.upper_bounds = current_iterate.feasibility_multipliers.upper_bounds + direction.feasibility_multipliers.upper_bounds;

      // there must be at least a nonzero dual to avoid trivial stationary points
      if (this->trial_multipliers.not_all_zero(model.number_variables, this->small_duals_threshold)) {
         // compute the ideal error (with a zero penalty parameter)
         const double infeasible_dual_error = l1Relaxation::compute_infeasible_dual_error(model, current_iterate);
         DEBUG << "Ideal dual error: " << infeasible_dual_error << '\n';
         const double scaled_error = infeasible_dual_error / std::max(1., current_iterate.primal_feasibility);
         this->penalty_parameter = std::min(this->penalty_parameter, scaled_error * scaled_error);
         DEBUG << "Further aggressively decrease the penalty parameter to " << this->penalty_parameter << '\n';
      }
      else {
         DEBUG << "l1Relaxation: all multipliers are almost 0. The penalty parameter won't be decreased\n";
      }
   }

   // measure that combines KKT error and complementarity error
   double l1Relaxation::compute_infeasible_dual_error(const Model& model, Iterate& current_iterate) const {
      const l1RelaxedProblem feasibility_problem{model, 0., this->constraint_violation_coefficient};
      // stationarity error
      feasibility_problem.evaluate_lagrangian_gradient(current_iterate.feasibility_residuals.lagrangian_gradient, current_iterate, this->trial_multipliers);
      double error = norm_1(current_iterate.residuals.lagrangian_gradient.constraints_contribution);

      // complementarity error
      constexpr double shift_value = 0.;
      error += feasibility_problem.complementarity_error(current_iterate.primals, current_iterate.evaluations.constraints,
            this->trial_multipliers, shift_value, Norm::L1);
      return error;
   }

   void l1Relaxation::enforce_linearized_residual_sufficient_decrease(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Direction& direction, double linearized_residual, double residual_lowest_violation, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      while (0. < this->penalty_parameter && !this->linearized_residual_sufficient_decrease(current_iterate, linearized_residual,
            residual_lowest_violation)) {
         // decrease the penalty parameter and re-solve the problem
         this->penalty_parameter /= this->parameters.decrease_factor;
         DEBUG << "Further decrease the penalty parameter to " << this->penalty_parameter << '\n';
         this->solve_l1_relaxed_problem(statistics, model, current_iterate, direction, this->penalty_parameter,
            trust_region_radius, warmstart_information);

         // recompute the linearized residual
         linearized_residual = model.constraint_violation(current_iterate.evaluations.constraints +
            current_iterate.evaluations.constraint_jacobian * direction.primals, Norm::L1);
         DEBUG << "Linearized infeasibility mk(dk): " << linearized_residual << "\n\n";
      }
      DEBUG << "Condition enforce_linearized_residual_sufficient_decrease is true\n";
   }

   bool l1Relaxation::linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual,
         double residual_lowest_violation) const {
      // if the feasibility problem is feasible wrt the original variables, look for a penalty parameter for which the l1 relaxed problem
      // also achieves feasibility
      if (residual_lowest_violation <= this->parameters.residual_small_threshold) {
         return (linearized_residual <= this->parameters.residual_small_threshold);
      }
      const double linearized_residual_reduction = current_iterate.progress.infeasibility - linearized_residual;
      const double lowest_linearized_residual_reduction = current_iterate.progress.infeasibility - residual_lowest_violation;
      if (lowest_linearized_residual_reduction < 0.) {
         WARNING << "lowest_linearized_residual_reduction = " << lowest_linearized_residual_reduction <<
         " is negative. Negative curvature in your problem\n";
      }
      return (linearized_residual_reduction >= this->parameters.epsilon1 * lowest_linearized_residual_reduction);
   }

   void l1Relaxation::enforce_descent_direction_for_l1_merit(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Direction& direction, const Direction& feasibility_direction, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      while (0. < this->penalty_parameter && !this->is_descent_direction_for_l1_merit_function(current_iterate, direction, feasibility_direction)) {
         // decrease the penalty parameter and re-solve the problem
         this->penalty_parameter /= this->parameters.decrease_factor;
         DEBUG << "Further decrease the penalty parameter to " << this->penalty_parameter << '\n';
         this->solve_l1_relaxed_problem(statistics, model, current_iterate, direction, this->penalty_parameter,
            trust_region_radius, warmstart_information);
      }
      DEBUG << "Condition enforce_descent_direction_for_l1_merit is true\n\n";
   }

   bool l1Relaxation::is_descent_direction_for_l1_merit_function(const Iterate& current_iterate, const Direction& direction,
         const Direction& feasibility_direction) const {
      const double predicted_l1_merit_reduction = current_iterate.primal_feasibility - direction.subproblem_objective;
      const double lowest_decrease_objective = current_iterate.primal_feasibility - feasibility_direction.subproblem_objective;
      return (predicted_l1_merit_reduction >= this->parameters.epsilon2 * lowest_decrease_objective);
   }

   bool l1Relaxation::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy, const Model& model,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      const l1RelaxedProblem l1_relaxed_problem{model, this->penalty_parameter, this->constraint_violation_coefficient};
      const bool accept_iterate = ConstraintRelaxationStrategy::is_iterate_acceptable(statistics, globalization_strategy,
         l1_relaxed_problem, *this->inequality_handling_method, current_iterate, trial_iterate, trial_iterate.multipliers,
         direction, step_length, user_callbacks);
      if (accept_iterate) {
         this->check_exact_relaxation(trial_iterate);
      }
      ConstraintRelaxationStrategy::set_primal_statistics(statistics, model, trial_iterate);
      warmstart_information.no_changes();
      return accept_iterate;
   }

   void l1Relaxation::compute_primal_dual_residuals(const Model& model, Iterate& iterate) {
      const l1RelaxedProblem l1_relaxed_problem{model, this->penalty_parameter, this->constraint_violation_coefficient};
      const l1RelaxedProblem feasibility_problem{model, 0., this->constraint_violation_coefficient};
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(model, l1_relaxed_problem, feasibility_problem, iterate);
   }

   void l1Relaxation::evaluate_progress_measures(InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, Iterate& iterate) const {
      this->set_infeasibility_measure(problem.model, iterate);
      this->set_objective_measure(problem.model, iterate);
      inequality_handling_method.set_auxiliary_measure(problem, iterate);
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

   std::string l1Relaxation::get_name() const {
      return "l1 relaxation " + this->inequality_handling_method->get_name() + " with " +
         this->l1_relaxed_hessian_model->get_name() + " Hessian and " + this->l1_relaxed_regularization_strategy->get_name() +
         " regularization";
   }

   size_t l1Relaxation::get_hessian_evaluation_count() const {
      return this->l1_relaxed_hessian_model->evaluation_count + this->feasibility_hessian_model->evaluation_count;
   }

   size_t l1Relaxation::get_number_subproblems_solved() const {
      return this->inequality_handling_method->number_subproblems_solved +
         this->feasibility_inequality_handling_method->number_subproblems_solved;
   }
} // namespace