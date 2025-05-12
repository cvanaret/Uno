// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/constraint_relaxation_strategies/l1RelaxedProblem.hpp"
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/regularization_strategies/UnstableRegularization.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   FeasibilityRestoration::FeasibilityRestoration(size_t number_bound_constraints, const Options& options) :
         ConstraintRelaxationStrategy(options),
         constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")),
         optimality_subproblem_layer(options),
         feasibility_subproblem_layer(options),
         optimality_inequality_handling_method(InequalityHandlingMethodFactory::create(number_bound_constraints, options)),
         feasibility_inequality_handling_method(InequalityHandlingMethodFactory::create(number_bound_constraints, options)),
         linear_feasibility_tolerance(options.get_double("tolerance")),
         switch_to_optimality_requires_linearized_feasibility(options.get_bool("switch_to_optimality_requires_linearized_feasibility")) {
   }

   void FeasibilityRestoration::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate,
         Direction& direction, const Options& options) {
      const OptimizationProblem optimality_problem{model};
      l1RelaxedProblem feasibility_problem{model, 0., this->constraint_violation_coefficient};
      this->reference_optimality_primals.resize(optimality_problem.number_variables);
      feasibility_problem.set_proximal_center(this->reference_optimality_primals.data());
      feasibility_problem.set_proximal_multiplier(this->feasibility_inequality_handling_method->proximal_coefficient());

      // memory allocation
      // TODO allocate the feasibility phase only when entering the first time?
      this->optimality_subproblem_layer.hessian_model->initialize(model);
      this->feasibility_subproblem_layer.hessian_model->initialize(model);
      this->optimality_inequality_handling_method->initialize(optimality_problem, *this->optimality_subproblem_layer.hessian_model,
         *this->optimality_subproblem_layer.regularization_strategy);
      this->feasibility_inequality_handling_method->initialize(feasibility_problem, *this->feasibility_subproblem_layer.hessian_model,
         *this->feasibility_subproblem_layer.regularization_strategy);
      direction = Direction(feasibility_problem.number_variables, feasibility_problem.number_constraints);

      // statistics
      this->optimality_subproblem_layer.initialize_statistics(statistics, options);
      this->feasibility_subproblem_layer.initialize_statistics(statistics, options);
      this->optimality_inequality_handling_method->initialize_statistics(statistics, options);
      this->feasibility_inequality_handling_method->initialize_statistics(statistics, options);
      statistics.add_column("phase", Statistics::int_width, options.get_int("statistics_restoration_phase_column_order"));
      statistics.set("phase", "OPT");

      // initial iterate
      initial_iterate.feasibility_residuals.lagrangian_gradient.resize(feasibility_problem.number_variables);
      initial_iterate.feasibility_multipliers.lower_bounds.resize(feasibility_problem.number_variables);
      initial_iterate.feasibility_multipliers.upper_bounds.resize(feasibility_problem.number_variables);
      this->optimality_inequality_handling_method->generate_initial_iterate(optimality_problem, initial_iterate);
      this->evaluate_progress_measures(*this->optimality_inequality_handling_method, model, initial_iterate);
      this->compute_primal_dual_residuals(model, initial_iterate);
      this->set_statistics(statistics, model, initial_iterate);
   }

   void FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, Direction& direction, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      direction.reset();
      // if we are in the optimality phase, solve the optimality problem
      if (this->current_phase == Phase::OPTIMALITY) {
         statistics.set("phase", "OPT");
         try {
            DEBUG << "Solving the optimality subproblem\n";
            const OptimizationProblem optimality_problem{model};
            this->solve_subproblem(statistics, *this->optimality_inequality_handling_method, optimality_problem, current_iterate,
               current_iterate.multipliers, direction, this->optimality_subproblem_layer, trust_region_radius, warmstart_information);
            if (direction.status == SubproblemStatus::INFEASIBLE) {
               // switch to the feasibility problem, starting from the current direction
               statistics.set("status", std::string("infeasible subproblem"));
               DEBUG << "/!\\ The subproblem is infeasible\n";
               this->switch_to_feasibility_problem(statistics, globalization_strategy, model, current_iterate, warmstart_information);
               this->feasibility_inequality_handling_method->set_initial_point(direction.primals);
            }
            else {
               warmstart_information.no_changes();
               return;
            }
         }
         catch (const UnstableRegularization&) {
            this->switch_to_feasibility_problem(statistics, globalization_strategy, model, current_iterate, warmstart_information);
         }
      }

      // solve the feasibility problem (minimize the constraint violation)
      DEBUG << "Solving the feasibility subproblem\n";
      statistics.set("phase", "FEAS");
      // note: failure of regularization should not happen here, since the feasibility Jacobian has full rank
      l1RelaxedProblem feasibility_problem{model, 0., this->constraint_violation_coefficient};
      feasibility_problem.set_proximal_center(this->reference_optimality_primals.data());
      feasibility_problem.set_proximal_multiplier(this->feasibility_inequality_handling_method->proximal_coefficient());
      this->solve_subproblem(statistics, *this->feasibility_inequality_handling_method, feasibility_problem, current_iterate,
         current_iterate.feasibility_multipliers, direction, this->feasibility_subproblem_layer, trust_region_radius,
         warmstart_information);
      std::swap(direction.multipliers, direction.feasibility_multipliers);
   }

   bool FeasibilityRestoration::solving_feasibility_problem() const {
      return (this->current_phase == Phase::FEASIBILITY_RESTORATION);
   }

   // precondition: this->current_phase == Phase::OPTIMALITY
   void FeasibilityRestoration::switch_to_feasibility_problem(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, WarmstartInformation& warmstart_information) {
      DEBUG << "Switching from optimality to restoration phase\n";
      this->current_phase = Phase::FEASIBILITY_RESTORATION;
      globalization_strategy.notify_switch_to_feasibility(current_iterate.progress);
      const l1RelaxedProblem feasibility_problem{model, 0., this->constraint_violation_coefficient};
      this->feasibility_inequality_handling_method->initialize_feasibility_problem(feasibility_problem, current_iterate);
      // save the current point (progress and primals) upon switching
      this->reference_optimality_progress = current_iterate.progress;
      this->reference_optimality_primals = current_iterate.primals;

      current_iterate.set_number_variables(feasibility_problem.number_variables);
      this->feasibility_inequality_handling_method->set_elastic_variable_values(feasibility_problem, current_iterate);
      DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

      if (Logger::level == INFO) statistics.print_current_line();
      warmstart_information.whole_problem_changed();
   }

   void FeasibilityRestoration::solve_subproblem(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers, Direction& direction,
         SubproblemLayer& subproblem_layer, double trust_region_radius, WarmstartInformation& warmstart_information) {
      direction.set_dimensions(problem.number_variables, problem.number_constraints);
      inequality_handling_method.solve(statistics, problem, current_iterate, current_multipliers, direction,
         subproblem_layer, trust_region_radius, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
   }

   bool FeasibilityRestoration::can_switch_to_optimality_phase(const Iterate& current_iterate, const GlobalizationStrategy& globalization_strategy,
         const Model& model, const Iterate& trial_iterate, const Direction& direction, double step_length) {
      return globalization_strategy.is_infeasibility_sufficiently_reduced(this->reference_optimality_progress, trial_iterate.progress) &&
         (!this->switch_to_optimality_requires_linearized_feasibility ||
         model.constraint_violation(current_iterate.evaluations.constraints + step_length*(current_iterate.evaluations.constraint_jacobian *
         direction.primals), this->residual_norm) <= this->linear_feasibility_tolerance);
   }

   void FeasibilityRestoration::switch_to_optimality_phase(Iterate& current_iterate, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& trial_iterate, WarmstartInformation& warmstart_information) {
      DEBUG << "Switching from restoration to optimality phase\n";
      this->current_phase = Phase::OPTIMALITY;
      globalization_strategy.notify_switch_to_optimality(current_iterate.progress);
      const OptimizationProblem optimality_problem{model};
      current_iterate.set_number_variables(optimality_problem.number_variables);
      trial_iterate.set_number_variables(optimality_problem.number_variables);
      current_iterate.objective_multiplier = trial_iterate.objective_multiplier = 1.;

      this->optimality_inequality_handling_method->exit_feasibility_problem(optimality_problem, trial_iterate);
      // set a cold start in the subproblem solver
      warmstart_information.whole_problem_changed();
   }

   bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      if (this->current_phase == Phase::OPTIMALITY) {
         const OptimizationProblem optimality_problem{model};
         this->optimality_inequality_handling_method->postprocess_iterate(optimality_problem, trial_iterate.primals, trial_iterate.multipliers);
         this->compute_progress_measures(*this->optimality_inequality_handling_method, model, globalization_strategy, current_iterate, trial_iterate);
      }
      else {
         const l1RelaxedProblem feasibility_problem{model, 0, this->constraint_violation_coefficient};
         this->feasibility_inequality_handling_method->postprocess_iterate(feasibility_problem, trial_iterate.primals, trial_iterate.feasibility_multipliers);
         this->compute_progress_measures(*this->feasibility_inequality_handling_method, model, globalization_strategy, current_iterate, trial_iterate);
      }
      trial_iterate.objective_multiplier = (this->current_phase == Phase::OPTIMALITY) ? 1. : 0.;
      const ProgressMeasures predicted_reduction = this->compute_predicted_reduction_models(
         (this->current_phase == Phase::OPTIMALITY) ? *this->optimality_inequality_handling_method : *this->feasibility_inequality_handling_method,
         model, current_iterate, direction, step_length);

      // possibly go from restoration phase to optimality phase
      if (this->current_phase == Phase::FEASIBILITY_RESTORATION && this->can_switch_to_optimality_phase(current_iterate, globalization_strategy,
            model, trial_iterate, direction, step_length)) {
         this->switch_to_optimality_phase(current_iterate, globalization_strategy, model, trial_iterate, warmstart_information);
      }
      else {
         warmstart_information.no_changes();
      }

      bool accept_iterate = false;
      const double objective_multiplier = (this->current_phase == Phase::OPTIMALITY) ? 1. : 0.;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_iterate.evaluate_objective(model);
         accept_iterate = true;
         statistics.set("status", "0 primal step");
      }
      else {
         accept_iterate = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reduction, objective_multiplier);
      }
      ConstraintRelaxationStrategy::set_progress_statistics(statistics, model, trial_iterate);
      if (accept_iterate) {
         user_callbacks.notify_acceptable_iterate(trial_iterate.primals,
               this->current_phase == Phase::OPTIMALITY ? trial_iterate.multipliers : trial_iterate.feasibility_multipliers,
               objective_multiplier);
      }
      return accept_iterate;
   }

   void FeasibilityRestoration::compute_primal_dual_residuals(const Model& model, Iterate& iterate) {
      const OptimizationProblem optimality_problem{model};
      const l1RelaxedProblem feasibility_problem{model, 0, this->constraint_violation_coefficient};
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(model, optimality_problem, feasibility_problem, iterate);
   }

   void FeasibilityRestoration::evaluate_progress_measures(InequalityHandlingMethod& inequality_handling_method, const Model& model, Iterate& iterate) const {
      this->set_infeasibility_measure(model, iterate);
      this->set_objective_measure(model, iterate);
      inequality_handling_method.set_auxiliary_measure(model, iterate);
   }

   ProgressMeasures FeasibilityRestoration::compute_predicted_reduction_models(InequalityHandlingMethod& inequality_handling_method,
         const Model& model, const Iterate& current_iterate, const Direction& direction, double step_length) const {
      return {
         this->compute_predicted_infeasibility_reduction(model, current_iterate, direction.primals, step_length),
         this->compute_predicted_objective_reduction(inequality_handling_method, current_iterate, direction.primals, step_length),
         inequality_handling_method.compute_predicted_auxiliary_reduction_model(model, current_iterate, direction.primals, step_length)
      };
   }

   void FeasibilityRestoration::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const {
      const auto& residuals = (this->current_phase == Phase::OPTIMALITY) ? iterate.residuals : iterate.feasibility_residuals;
      statistics.set("stationarity", residuals.stationarity);
      statistics.set("complementarity", residuals.complementarity);
   }

   std::string FeasibilityRestoration::get_name() const {
      return "feasibility restoration";
   }

   size_t FeasibilityRestoration::get_hessian_evaluation_count() const {
      return this->optimality_subproblem_layer.get_hessian_evaluation_count() +
         this->feasibility_subproblem_layer.get_hessian_evaluation_count();
   }

   size_t FeasibilityRestoration::get_number_subproblems_solved() const {
      return this->optimality_inequality_handling_method->number_subproblems_solved +
         this->feasibility_inequality_handling_method->number_subproblems_solved;
   }
} // namespace