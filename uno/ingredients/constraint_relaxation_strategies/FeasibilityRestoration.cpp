// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <functional>
#include "FeasibilityRestoration.hpp"
#include "relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategyFactory.hpp"
#include "linear_algebra/VectorView.hpp"
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

namespace uno {
   FeasibilityRestoration::FeasibilityRestoration(const Model& model, bool use_trust_region, const Options& options) :
         ConstraintRelaxationStrategy(options),
         constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")),
         original_problem(model),
         // relax the linear constraints in the l1 relaxed problem only if we are using a trust-region constraint
         feasibility_problem(model, 0., this->constraint_violation_coefficient, use_trust_region),
         hessian_model(HessianModelFactory::create(model, 1., options)),
         feasibility_hessian_model(HessianModelFactory::create(model, 0., options)),
         inertia_correction_strategy(InertiaCorrectionStrategyFactory::create(options)),
         feasibility_inertia_correction_strategy(InertiaCorrectionStrategyFactory::create(options)),
         inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         feasibility_inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         globalization_strategy(GlobalizationStrategyFactory::create(model, options)),
         feasibility_globalization_strategy(options),
         linear_feasibility_tolerance(options.get_double("primal_tolerance")),
         switch_to_optimality_requires_linearized_feasibility(options.get_bool("switch_to_optimality_requires_linearized_feasibility")) {
   }

   void FeasibilityRestoration::initialize(Statistics& statistics, Iterate& initial_iterate, Direction& direction,
         bool uses_trust_region, EvaluationCache& evaluation_cache) {
      this->reference_optimality_primals.resize(this->original_problem.number_variables);

      // memory allocation
      this->inequality_handling_method->initialize(this->original_problem, initial_iterate,
         *this->hessian_model, *this->inertia_correction_strategy, uses_trust_region);
      direction = Direction(
         std::max(this->original_problem.number_variables, this->feasibility_problem.number_variables),
         std::max(this->original_problem.number_constraints, this->feasibility_problem.number_constraints)
      );

      // statistics
      this->inertia_correction_strategy->initialize_statistics(statistics);
      this->inequality_handling_method->initialize_statistics(statistics);
      this->feasibility_inertia_correction_strategy->initialize_statistics(statistics);
      this->feasibility_inequality_handling_method->initialize_statistics(statistics);
      this->hessian_model->initialize_statistics(statistics);
      this->feasibility_hessian_model->initialize_statistics(statistics);
      statistics.add_column("Phase", Statistics::int_width - 1, 3, Statistics::column_order.at("Phase"));
      statistics.set("Phase", "OPT");

      // initial iterate
      this->inequality_handling_method->generate_initial_iterate(initial_iterate, evaluation_cache);
      this->compute_residuals(this->original_problem, initial_iterate, evaluation_cache.current_evaluations);
      this->globalization_strategy->initialize(statistics, initial_iterate);
      this->feasibility_globalization_strategy.initialize(statistics, initial_iterate);
   }

   void FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double trust_region_radius, Evaluations& current_evaluations, WarmstartInformation& warmstart_information) {
      direction.reset();
      const bool use_trust_region = (trust_region_radius < INF<double>);

      // if we are in the optimality phase, solve the optimality problem
      if (this->current_phase == Phase::OPTIMALITY) {
         DEBUG << "Solving the optimality subproblem\n";
         statistics.set("Phase", "OPT");
         this->solve_subproblem(statistics, *this->inequality_handling_method, this->original_problem,
            current_iterate, direction, trust_region_radius, current_evaluations, warmstart_information);
         if (direction.status == SubproblemStatus::INFEASIBLE) {
            // switch to the feasibility problem, starting from the current direction
            statistics.set("Status", std::string("infeasible"));
            DEBUG << "/!\\ The subproblem is infeasible\n";
            this->switch_to_feasibility_problem(statistics, current_iterate, current_evaluations, use_trust_region,
               warmstart_information);
            this->feasibility_inequality_handling_method->set_initial_point(direction.primals);
         }
         else {
            warmstart_information.no_changes();
            return;
         }
      }

      // solve the feasibility problem (minimize the constraint violation)
      DEBUG << "Solving the feasibility subproblem\n";
      statistics.set("Phase", "FEAS");
      // note: failure of regularization should not happen here, since the feasibility Jacobian has full rank
      this->feasibility_problem.set_proximal_coefficient(this->inequality_handling_method->proximal_coefficient());
      this->solve_subproblem(statistics, *this->feasibility_inequality_handling_method, this->feasibility_problem,
         current_iterate, direction, trust_region_radius, current_evaluations, warmstart_information);
   }

   bool FeasibilityRestoration::solving_feasibility_problem() const {
      return (this->current_phase == Phase::FEASIBILITY_RESTORATION);
   }

   // precondition: this->current_phase == Phase::OPTIMALITY
   void FeasibilityRestoration::switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         Evaluations& current_evaluations, bool uses_trust_region, WarmstartInformation& warmstart_information) {
      DEBUG << "\nSwitching from optimality to restoration phase\n";
      this->current_phase = Phase::FEASIBILITY_RESTORATION;
      this->globalization_strategy->notify_switch_to_feasibility(current_iterate.progress);

      // save the current point (progress and primals) upon switching
      this->reference_optimality_progress = current_iterate.progress;
      this->reference_optimality_primals = current_iterate.primals;
      this->feasibility_problem.set_proximal_coefficient(this->inequality_handling_method->proximal_coefficient());
      this->feasibility_problem.set_proximal_center(this->reference_optimality_primals.data());

      current_iterate.set_number_variables(this->feasibility_problem.number_variables);
      // swap the iterate's multipliers and the feasibility multipliers maintained by the class
      this->other_phase_multipliers.constraints.resize(this->feasibility_problem.number_constraints);
      this->other_phase_multipliers.lower_bounds.resize(this->feasibility_problem.number_variables);
      this->other_phase_multipliers.upper_bounds.resize(this->feasibility_problem.number_variables);
      std::swap(current_iterate.multipliers, this->other_phase_multipliers);

      this->feasibility_inequality_handling_method->initialize_feasibility_problem(current_iterate);
      this->feasibility_inequality_handling_method->set_elastic_variable_values(this->feasibility_problem, current_iterate,
         current_evaluations);

      DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

      // initialize the feasibility ingredients upon the first switch to feasibility restoration
      if (this->first_switch_to_feasibility) {
         this->feasibility_inequality_handling_method->initialize(this->feasibility_problem, current_iterate,
            *this->feasibility_hessian_model, *this->feasibility_inertia_correction_strategy, uses_trust_region);
         this->first_switch_to_feasibility = false;
      }

      if (Logger::level == INFO) statistics.print_current_line();
      warmstart_information.whole_problem_changed();
   }

   void FeasibilityRestoration::solve_subproblem(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, Iterate& current_iterate, Direction& direction, double trust_region_radius,
         Evaluations& current_evaluations, WarmstartInformation& warmstart_information) {
      direction.set_dimensions(problem.number_variables, problem.number_constraints);
      inequality_handling_method.solve(statistics, current_iterate, direction, trust_region_radius, current_evaluations,
         warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
   }

   bool FeasibilityRestoration::can_switch_to_optimality_phase(const Model& model, const Iterate& trial_iterate,
         const Direction& direction, double step_length, EvaluationCache& evaluation_cache) const {
      if (this->globalization_strategy->is_infeasibility_sufficiently_reduced(this->reference_optimality_progress,
            trial_iterate.progress)) {
         if (!this->switch_to_optimality_requires_linearized_feasibility) {
            return true;
         }
         // compute the linearized constraint violation
         // TODO preallocate
         Vector<double> result(model.number_constraints);
         evaluation_cache.current_evaluations.compute_jacobian_vector_product(direction.primals, result);
         const double trial_linearized_constraint_violation = model.constraint_violation(evaluation_cache.current_evaluations.constraints +
            step_length * result, this->residual_norm);
         return (trial_linearized_constraint_violation <= this->linear_feasibility_tolerance);
      }
      return false;
   }

   void FeasibilityRestoration::switch_back_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate) {
      DEBUG << "Switching from restoration back to optimality phase\n";
      this->current_phase = Phase::OPTIMALITY;
      this->globalization_strategy->notify_switch_to_optimality(current_iterate.progress);

      current_iterate.set_number_variables(this->original_problem.number_variables);
      // swap the iterate's multipliers and the optimality multipliers maintained by the class
      std::swap(current_iterate.multipliers, this->other_phase_multipliers);
      trial_iterate.set_number_variables(this->original_problem.number_variables);
      current_iterate.objective_multiplier = trial_iterate.objective_multiplier = 1.;
   }

   bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, const Model& model,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         EvaluationCache& evaluation_cache, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      bool accept_iterate = false;
      // determine acceptability, depending on the current phase
      if (this->current_phase == Phase::OPTIMALITY) {
         accept_iterate = this->inequality_handling_method->is_iterate_acceptable(statistics, *this->globalization_strategy,
            current_iterate, trial_iterate, direction, step_length, evaluation_cache, user_callbacks);
         if (accept_iterate) {
            this->hessian_model->notify_accepted_iterate(statistics, current_iterate, trial_iterate, evaluation_cache);
         }
      }
      else {
         accept_iterate = this->feasibility_inequality_handling_method->is_iterate_acceptable(statistics,
            this->feasibility_globalization_strategy, current_iterate, trial_iterate, direction, step_length,
            evaluation_cache, user_callbacks);
         if (accept_iterate) {
            this->feasibility_hessian_model->notify_accepted_iterate(statistics, current_iterate, trial_iterate, evaluation_cache);
         }
      }

      // possibly go from restoration phase to optimality phase
      if (trial_iterate.status == SolutionStatus::NOT_OPTIMAL && this->current_phase == Phase::FEASIBILITY_RESTORATION &&
            this->can_switch_to_optimality_phase(model, trial_iterate, direction, step_length, evaluation_cache)) {
         this->switch_back_to_optimality_phase(current_iterate, trial_iterate);
         // set a cold start in the subproblem solver
         warmstart_information.whole_problem_changed();
      }
      else {
         warmstart_information.no_changes();
      }

      // check termination
      if (this->current_phase == Phase::OPTIMALITY) {
         this->compute_residuals(this->original_problem, trial_iterate, evaluation_cache.trial_evaluations);
         trial_iterate.status = this->check_termination(this->original_problem, trial_iterate,
            evaluation_cache.trial_evaluations);
      }
      else {
         this->compute_residuals(this->feasibility_problem, trial_iterate, evaluation_cache.trial_evaluations);
         trial_iterate.status = this->check_termination(this->feasibility_problem, trial_iterate,
            evaluation_cache.trial_evaluations);
      }
      return accept_iterate;
   }

   std::string FeasibilityRestoration::get_name() const {
      return this->globalization_strategy->get_name() + " restoration " + this->inequality_handling_method->get_name() +
         " with " + this->hessian_model->name + " Hessian and " + this->inertia_correction_strategy->get_name() +
         " regularization";
   }

   size_t FeasibilityRestoration::get_number_subproblems_solved() const {
      return this->inequality_handling_method->number_subproblems_solved +
         this->feasibility_inequality_handling_method->number_subproblems_solved;
   }
} // namespace