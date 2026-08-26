// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <functional>
#include "FeasibilityRestoration.hpp"
#include "relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "linear_algebra/View.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Sum.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   FeasibilityRestoration::FeasibilityRestoration(const Model& model, bool use_trust_region, Options& options) :
         ConstraintRelaxationStrategy(options),
         constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")),
         original_problem(model),
         // relax the linear constraints in the l1 relaxed problem only if we are using a trust-region constraint
         feasibility_problem(model, 0., this->constraint_violation_coefficient, use_trust_region),
         globalization_strategy(GlobalizationStrategyFactory::create(model, options)),
         feasibility_globalization_strategy(options),
         linear_feasibility_tolerance(options.get_double("primal_tolerance")),
         switch_to_optimality_requires_linearized_feasibility(options.get_bool("switch_to_optimality_requires_linearized_feasibility")) {
   }

   FeasibilityRestoration::~FeasibilityRestoration() = default;

   void FeasibilityRestoration::initialize(Statistics& statistics, Iterate& initial_iterate, bool uses_trust_region,
         EvaluationCache& evaluation_cache, Options& options) {
      this->initial_point.resize(this->original_problem.number_variables);
      this->reference_optimality_primals.resize(this->original_problem.number_variables);
      this->feasibility_problem.set_proximal_center(this->reference_optimality_primals.data());

      // reformulation of the original problem and the feasibility problem
      INFO << "- Allocating optimality (original) problem:\n";
      this->inequality_handling_method = InequalityHandlingMethodFactory::create(this->original_problem, uses_trust_region,
         1., options);
      INFO << "- Allocating feasibility problem:\n";
      this->feasibility_inequality_handling_method = InequalityHandlingMethodFactory::create(this->feasibility_problem,
         uses_trust_region, 0., options);

      // initial iterate
      this->inequality_handling_method->generate_initial_iterate(initial_iterate, evaluation_cache.current_evaluations);
      this->inequality_handling_method->evaluate_progress_measures(initial_iterate, evaluation_cache.current_evaluations);
      this->compute_residuals(this->original_problem, initial_iterate, evaluation_cache.current_evaluations);
      this->globalization_strategy->initialize(statistics, initial_iterate);
      this->feasibility_globalization_strategy.initialize(statistics, initial_iterate);

      // statistics
      this->inequality_handling_method->initialize_statistics(statistics);
      this->feasibility_inequality_handling_method->initialize_statistics(statistics);
      statistics.add_column("Phase", Statistics::int_width - 1, 3);
      statistics.set("Phase", "OPT");
   }

   const Direction& FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
         double trust_region_radius, Evaluations& current_evaluations, WarmstartInformation& warmstart_information) {

      // if we are in the optimality phase, solve the optimality problem
      if (this->current_phase == Phase::OPTIMALITY) {
         DEBUG << "Solving the optimality subproblem\n";
         statistics.set("Phase", "OPT");
         this->initial_point.fill(0.);
         const Direction& optimality_direction = this->solve_subproblem(statistics, *this->inequality_handling_method,
            *this->globalization_strategy, current_iterate, trust_region_radius, current_evaluations, warmstart_information);
         if (optimality_direction.status == SubproblemStatus::INFEASIBLE) {
            // switch to the feasibility problem, starting from the current direction
            statistics.set("Status", std::string("infeasible"));
            DEBUG << "/!\\ The subproblem is infeasible\n";
            this->initial_point = view(optimality_direction.primals, 0, this->original_problem.number_variables);
            this->switch_to_feasibility_problem(statistics, current_iterate, current_evaluations, warmstart_information);
         }
         else {
            warmstart_information.no_changes();
            return optimality_direction;
         }
      }

      // solve the feasibility problem (minimize the constraint violation)
      DEBUG << "Solving the feasibility subproblem\n";
      statistics.set("Phase", "FEAS");
      // note: failure of regularization should not happen here, since the feasibility Jacobian has full rank
      const Direction& feasibility_direction = this->solve_subproblem(statistics, *this->feasibility_inequality_handling_method,
         this->feasibility_globalization_strategy, current_iterate, trust_region_radius, current_evaluations,
         warmstart_information);
      return feasibility_direction;
   }

   bool FeasibilityRestoration::solving_feasibility_problem() const {
      return (this->current_phase == Phase::FEASIBILITY_RESTORATION);
   }

   // precondition: this->current_phase == Phase::OPTIMALITY
   void FeasibilityRestoration::switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         Evaluations& current_evaluations, WarmstartInformation& warmstart_information) {
      DEBUG << "\nSwitching from optimality to restoration phase\n";
      this->current_phase = Phase::FEASIBILITY_RESTORATION;
      this->globalization_strategy->notify_switch_to_feasibility(current_iterate.progress);

      // save the current point (progress and primals) upon switching
      this->reference_optimality_progress = current_iterate.progress;
      this->reference_optimality_primals = view(current_iterate.primals, 0, this->original_problem.number_variables);
      this->feasibility_problem.set_proximal_coefficient(this->inequality_handling_method->proximal_coefficient());
      this->feasibility_problem.set_proximal_center(this->reference_optimality_primals.data());

      this->initial_point.resize(this->feasibility_problem.number_variables);
      // swap the iterate's multipliers and the feasibility multipliers maintained by the class
      if (this->first_switch_to_feasibility) {
         this->other_phase_multipliers.constraints.resize(this->feasibility_problem.number_constraints);
         this->other_phase_multipliers.lower_bounds.resize(this->feasibility_problem.number_variables);
         this->other_phase_multipliers.upper_bounds.resize(this->feasibility_problem.number_variables);
         this->first_switch_to_feasibility = false;
      }
      std::swap(current_iterate.multipliers, this->other_phase_multipliers);

      const auto [number_variables_optimality, _] = this->inequality_handling_method->get_problem_dimensions();
      const auto [number_variables_feasibility, __] = this->feasibility_inequality_handling_method->get_problem_dimensions();

      const auto auxiliary_variables = view(current_iterate.primals, this->original_problem.number_variables,
         number_variables_optimality);
      const auto auxiliary_variables_lb_duals = view(current_iterate.multipliers.lower_bounds, this->original_problem.number_variables,
         number_variables_optimality);
      const auto auxiliary_variables_ub_duals = view(current_iterate.multipliers.lower_bounds, this->original_problem.number_variables,
         number_variables_optimality);
      this->feasibility_inequality_handling_method->initialize_feasibility_problem(current_iterate);
      // copy the additional variables into the last block
      view(current_iterate.primals, this->feasibility_problem.number_variables, number_variables_feasibility) = auxiliary_variables;
      view(current_iterate.multipliers.lower_bounds, this->feasibility_problem.number_variables, number_variables_feasibility) = auxiliary_variables_lb_duals;
      view(current_iterate.multipliers.upper_bounds, this->feasibility_problem.number_variables, number_variables_feasibility) = auxiliary_variables_ub_duals;
      this->feasibility_inequality_handling_method->set_elastic_variable_values(this->feasibility_problem, current_iterate,
         current_evaluations);
      // re-evaluate the progress measures at the current iterate
      this->feasibility_inequality_handling_method->evaluate_progress_measures(current_iterate, current_evaluations);

      DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

      if (Logger::level == INFO) statistics.print_current_line();
      warmstart_information.whole_problem_changed();
   }

   bool FeasibilityRestoration::has_second_order_corrections() const {
      if (this->current_phase == Phase::OPTIMALITY) {
         return this->inequality_handling_method->has_second_order_corrections();
      }
      else {
         return this->feasibility_inequality_handling_method->has_second_order_corrections();
      }
   }

   void FeasibilityRestoration::initialize_second_order_corrections(const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      if (this->current_phase == Phase::OPTIMALITY) {
         this->inequality_handling_method->initialize_second_order_corrections(current_iterate, trial_iterate, current_evaluations,
            trial_evaluations);
      }
      else {
         this->feasibility_inequality_handling_method->initialize_second_order_corrections(current_iterate, trial_iterate,
            current_evaluations, trial_evaluations);
      }
   }

   const Direction& FeasibilityRestoration::compute_second_order_correction(Iterate& current_iterate) {
      if (this->current_phase == Phase::OPTIMALITY) {
         return this->inequality_handling_method->compute_second_order_correction(current_iterate);
      }
      else {
         return this->feasibility_inequality_handling_method->compute_second_order_correction(current_iterate);
      }
   }

   void FeasibilityRestoration::update_second_order_corrections(const Iterate& trial_iterate, Evaluations& trial_evaluations) {
      if (this->current_phase == Phase::OPTIMALITY) {
         return this->inequality_handling_method->update_second_order_corrections(trial_iterate, trial_evaluations);
      }
      else {
         return this->feasibility_inequality_handling_method->update_second_order_corrections(trial_iterate, trial_evaluations);
      }
   }

   const Direction& FeasibilityRestoration::solve_subproblem(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method,
         GlobalizationStrategy& globalization_strategy, Iterate& current_iterate, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      // if the problem definition changed, reset the globalization strategy and recompute the current auxiliary measure
      if (inequality_handling_method.update_parameterization(statistics, current_iterate)) {
         globalization_strategy.reset();
         inequality_handling_method.evaluate_progress_measures(current_iterate, current_evaluations); // TODO auxiliary
      }

      const Direction& direction = inequality_handling_method.solve(statistics, current_iterate, trust_region_radius,
         this->initial_point, current_evaluations, warmstart_information);
      ++this->number_subproblems_solved;
      this->initial_point.fill(0.);
      DEBUG3 << direction << '\n';
      return direction;
   }

   bool FeasibilityRestoration::can_switch_to_optimality_phase(const Model& model, const Iterate& trial_iterate,
         const Direction& direction, double step_length, Evaluations& current_evaluations) const {
      if (this->globalization_strategy->is_infeasibility_sufficiently_reduced(this->reference_optimality_progress,
            trial_iterate.progress)) {
         if (!this->switch_to_optimality_requires_linearized_feasibility) {
            return true;
         }
         // compute the linearized constraint violation
         // TODO preallocate
         Vector<double> result(model.number_constraints);
         current_evaluations.compute_jacobian_vector_product(model, direction.primals.data(), result.data());
         const double trial_linearized_constraint_violation = model.constraint_violation(current_evaluations.constraints +
            step_length * result, this->residual_norm);
         return (trial_linearized_constraint_violation <= this->linear_feasibility_tolerance);
      }
      return false;
   }

   void FeasibilityRestoration::switch_back_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      DEBUG << "Switching from restoration back to optimality phase\n";
      this->current_phase = Phase::OPTIMALITY;
      this->inequality_handling_method->evaluate_progress_measures(current_iterate, current_evaluations);
      this->inequality_handling_method->evaluate_progress_measures(trial_iterate, trial_evaluations);
      this->globalization_strategy->notify_switch_to_optimality(current_iterate.progress);

      // swap the iterate's multipliers and the optimality multipliers maintained by the class, and possibly compute
      // least-squares multipliers for the original problem
      std::swap(current_iterate.multipliers, this->other_phase_multipliers);
      this->inequality_handling_method->compute_least_squares_multipliers(trial_iterate, trial_evaluations);

      const auto [number_variables_optimality, _] = this->inequality_handling_method->get_problem_dimensions();
      const auto [number_variables_feasibility, __] = this->feasibility_inequality_handling_method->get_problem_dimensions();

      const auto current_auxiliary_variables = view(current_iterate.primals, this->feasibility_problem.number_variables,
         number_variables_feasibility);
      const auto current_auxiliary_variables_lb_duals = view(current_iterate.multipliers.lower_bounds, this->feasibility_problem.number_variables,
         number_variables_feasibility);
      const auto current_auxiliary_variables_ub_duals = view(current_iterate.multipliers.upper_bounds, this->feasibility_problem.number_variables,
         number_variables_feasibility);
      view(current_iterate.primals, this->original_problem.number_variables, number_variables_optimality) = current_auxiliary_variables;
      view(current_iterate.multipliers.lower_bounds, this->original_problem.number_variables, number_variables_optimality) = current_auxiliary_variables_lb_duals;
      view(current_iterate.multipliers.upper_bounds, this->original_problem.number_variables, number_variables_optimality) = current_auxiliary_variables_ub_duals;
      current_iterate.set_number_variables(number_variables_optimality);

      const auto trial_auxiliary_variables = view(trial_iterate.primals, this->feasibility_problem.number_variables,
         number_variables_feasibility);
      const auto trial_auxiliary_variables_lb_duals = view(trial_iterate.multipliers.lower_bounds, this->feasibility_problem.number_variables,
         number_variables_feasibility);
      const auto trial_auxiliary_variables_ub_duals = view(trial_iterate.multipliers.upper_bounds, this->feasibility_problem.number_variables,
         number_variables_feasibility);
      view(trial_iterate.primals, this->original_problem.number_variables, number_variables_optimality) = trial_auxiliary_variables;
      view(trial_iterate.multipliers.lower_bounds, this->original_problem.number_variables, number_variables_optimality) = trial_auxiliary_variables_lb_duals;
      view(trial_iterate.multipliers.upper_bounds, this->original_problem.number_variables, number_variables_optimality) = trial_auxiliary_variables_ub_duals;
      trial_iterate.set_number_variables(number_variables_optimality);

      current_iterate.objective_multiplier = trial_iterate.objective_multiplier = 1.;
      this->initial_point.resize(this->original_problem.number_variables);
   }

   bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, const Model& model,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         bool uses_trust_region, Evaluations& current_evaluations, Evaluations& trial_evaluations,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      bool accept_iterate = false;
      // determine acceptability, depending on the current phase
      if (this->current_phase == Phase::OPTIMALITY) {
         accept_iterate = this->inequality_handling_method->is_iterate_acceptable(statistics, *this->globalization_strategy,
            current_iterate, trial_iterate, direction, step_length, current_evaluations, trial_evaluations);
         if (uses_trust_region || accept_iterate) {
            this->inequality_handling_method->notify_trial_iterate(statistics, current_iterate, trial_iterate, current_evaluations,
               trial_evaluations);
         }
      }
      else {
         accept_iterate = this->feasibility_inequality_handling_method->is_iterate_acceptable(statistics,
            this->feasibility_globalization_strategy, current_iterate, trial_iterate, direction, step_length, current_evaluations,
            trial_evaluations);
         if (uses_trust_region || accept_iterate) {
            this->feasibility_inequality_handling_method->notify_trial_iterate(statistics, current_iterate, trial_iterate,
               current_evaluations, trial_evaluations);
         }
      }

      // possibly go from restoration phase to optimality phase
      if (accept_iterate && this->current_phase == Phase::FEASIBILITY_RESTORATION && this->can_switch_to_optimality_phase(model,
            trial_iterate, direction, step_length, current_evaluations)) {
         this->switch_back_to_optimality_phase(current_iterate, trial_iterate, current_evaluations, trial_evaluations);
         // set a cold start in the subproblem solver
         warmstart_information.whole_problem_changed();
      }
      else {
         warmstart_information.no_changes();
      }

      // check termination
      if (this->current_phase == Phase::OPTIMALITY) {
         this->compute_residuals(this->original_problem, trial_iterate, trial_evaluations);
         trial_iterate.status = this->check_termination(this->original_problem, trial_iterate, trial_evaluations);
         if (accept_iterate) {
            user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers,
               this->original_problem.get_objective_multiplier(), trial_iterate.progress.infeasibility,
               trial_iterate.residuals.stationarity, trial_iterate.residuals.complementarity);
         }
      }
      else {
         this->compute_residuals(this->feasibility_problem, trial_iterate, trial_evaluations);
         trial_iterate.status = this->check_termination(this->feasibility_problem, trial_iterate, trial_evaluations);
         if (accept_iterate) {
            user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers,
               this->feasibility_problem.get_objective_multiplier(), trial_iterate.progress.infeasibility,
               trial_iterate.residuals.stationarity, trial_iterate.residuals.complementarity);
         }
      }
      return accept_iterate;
   }

   std::string FeasibilityRestoration::get_name() const {
      return this->globalization_strategy->get_name() + " restoration " + this->inequality_handling_method->get_name();
   }
} // namespace