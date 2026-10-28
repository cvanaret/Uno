// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <memory>
#include "NoRelaxation.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   NoRelaxation::NoRelaxation(const Model& model, const Options& options):
         ConstraintRelaxationStrategy(options),
         original_problem(model, options.get_double("residual_scaling_threshold")),
         globalization_strategy(options),
         least_square_multiplier_max_norm(options.get_double("least_square_multiplier_max_norm")) {
   }

   NoRelaxation::~NoRelaxation() = default;

   void NoRelaxation::initialize(Statistics& statistics, Iterate& initial_iterate, bool uses_trust_region,
         EvaluationCache& evaluation_cache, const Options& options, std::vector<OptionOverride>& option_overrides) {
      this->initial_point.resize(this->original_problem.number_variables);

      // reformulation of the original problem
      INFO << "- Allocating method: ";
      this->inequality_handling_method = InequalityHandlingMethodFactory::create(this->original_problem, uses_trust_region,
         1., options, option_overrides);
      this->inequality_handling_method->initialize_memory();
      INFO << '\n';

      // initial iterate
      this->inequality_handling_method->create_initial_iterate(initial_iterate, evaluation_cache.current_evaluations,
         this->least_square_multiplier_max_norm);
      this->inequality_handling_method->evaluate_progress_measures(initial_iterate, evaluation_cache.current_evaluations);
      this->inequality_handling_method->compute_residuals(initial_iterate, evaluation_cache.current_evaluations);
      this->globalization_strategy.initialize(statistics, initial_iterate);

      // statistics
      this->inequality_handling_method->initialize_statistics(statistics);
   }

   void NoRelaxation::postprocess_solution(Iterate& iterate, Evaluations& evaluations) const {
      this->inequality_handling_method->postprocess_solution(iterate, evaluations);
   }

   const Direction& NoRelaxation::compute_direction(Statistics& statistics, Iterate& current_iterate,
         double trust_region_radius, Evaluations& current_evaluations, WarmstartInformation& warmstart_information) {
      DEBUG << "Solving the subproblem\n";
      const bool parameterization_updated = this->inequality_handling_method->update_parameterization(statistics,
         current_iterate, current_evaluations);
      // if the problem definition changed, reset the globalization strategy and recompute the current auxiliary measure
      if (parameterization_updated) {
         ++this->number_parameterization_updates;
         this->globalization_strategy.reset();
         this->inequality_handling_method->evaluate_progress_measures(current_iterate, current_evaluations); // TODO auxiliary
      }

      this->initial_point.fill(0.);
      const Direction& direction = this->inequality_handling_method->solve(statistics, current_iterate, trust_region_radius,
         this->initial_point, current_evaluations, warmstart_information);
      ++this->number_subproblems_solved;
      DEBUG3 << direction << '\n';
      warmstart_information.no_changes();
      return direction;
   }

   bool NoRelaxation::solving_feasibility_problem() const {
      return false;
   }

   void NoRelaxation::switch_to_feasibility_problem(Statistics& /*statistics*/, Iterate& /*current_iterate*/,
         Evaluations& /*current_evaluations*/, WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("Switching to the feasibility problem should not happen");
   }

   bool NoRelaxation::has_second_order_corrections() const {
      return this->inequality_handling_method->has_second_order_corrections();
   }

   void NoRelaxation::initialize_second_order_corrections(const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      this->inequality_handling_method->initialize_second_order_corrections(current_iterate, trial_iterate, current_evaluations,
         trial_evaluations);
   }

   const Direction& NoRelaxation::compute_second_order_correction(Statistics& statistics, Iterate& current_iterate) {
      return this->inequality_handling_method->compute_second_order_correction(statistics, current_iterate);
   }

   void NoRelaxation::update_second_order_corrections(const Iterate& trial_iterate, Evaluations& trial_evaluations) {
      this->inequality_handling_method->update_second_order_corrections(trial_iterate, trial_evaluations);
   }

   PredictedReductionModels NoRelaxation::build_predicted_reduction_models(const Iterate& current_iterate,
         const Direction& direction, Evaluations& current_evaluations) const {
      return this->inequality_handling_method->build_predicted_reduction_models(current_iterate, direction, current_evaluations);
   }

   Decision NoRelaxation::is_iterate_acceptable(Statistics& statistics, const Model& /*model*/, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double /*step_length*/, bool uses_trust_region, Evaluations& current_evaluations,
         Evaluations& trial_evaluations, const ProgressMeasures& predicted_reductions, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) {
      const Decision decision = this->inequality_handling_method->is_iterate_acceptable(statistics, this->globalization_strategy,
         current_iterate, trial_iterate, direction, trial_evaluations, predicted_reductions);
      if (decision.accepted) {
         this->inequality_handling_method->compute_residuals(trial_iterate, trial_evaluations);
         trial_iterate.status = this->check_termination(this->original_problem, trial_iterate, trial_evaluations);
         user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers,
            this->original_problem.get_objective_multiplier(), trial_iterate.scaled_residuals.primal_infeasibility,
            trial_iterate.scaled_residuals.stationarity, trial_iterate.scaled_residuals.complementarity);
      }
      if (uses_trust_region || decision.accepted) {
         this->inequality_handling_method->notify_trial_iterate(statistics, current_iterate, trial_iterate, current_evaluations,
            trial_evaluations);
      }
      warmstart_information.no_changes();
      return decision;
   }

   void NoRelaxation::notify_step_limited_by(RejectionReason sole_rejection_reason) {
      this->globalization_strategy.notify_step_limited_by(sole_rejection_reason);
   }

   bool NoRelaxation::accept_iterate_unconditionally(Statistics& statistics, const Model& /*model*/, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& /*direction*/, double /*step_length*/, bool force_parameterization_update,
         Evaluations& current_evaluations, Evaluations& trial_evaluations, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) {
      // update the parameterization first: the progress measures and the postprocessing of the trial iterate depend on it
      const bool parameterization_updated = force_parameterization_update &&
         this->inequality_handling_method->force_parameterization_update(statistics);
      if (parameterization_updated) {
         ++this->number_parameterization_updates;
         this->globalization_strategy.reset();
         this->inequality_handling_method->evaluate_progress_measures(current_iterate, current_evaluations);
      }
      this->inequality_handling_method->accept_iterate_unconditionally(trial_iterate, trial_evaluations);
      this->inequality_handling_method->compute_residuals(trial_iterate, trial_evaluations);
      trial_iterate.status = this->check_termination(this->original_problem, trial_iterate, trial_evaluations);
      user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers,
         this->original_problem.get_objective_multiplier(), trial_iterate.scaled_residuals.primal_infeasibility,
         trial_iterate.scaled_residuals.stationarity, trial_iterate.scaled_residuals.complementarity);
      this->inequality_handling_method->notify_trial_iterate(statistics, current_iterate, trial_iterate, current_evaluations,
         trial_evaluations);
      warmstart_information.no_changes();
      return parameterization_updated;
   }

   bool NoRelaxation::terminate_at_line_search_failure(const Iterate& current_iterate, Iterate& trial_iterate,
         const Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      return this->fall_back_after_line_search_failure(this->original_problem, current_iterate, trial_iterate,
         current_evaluations, trial_evaluations);
   }

   double NoRelaxation::small_infeasibility_threshold() const {
      return this->globalization_strategy.small_infeasibility_threshold();
   }

   std::string NoRelaxation::get_name() const {
      return this->globalization_strategy.get_name() + " " + this->inequality_handling_method->get_name();
   }
} // namespace