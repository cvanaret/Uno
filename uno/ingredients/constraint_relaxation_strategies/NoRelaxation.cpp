// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <memory>
#include "NoRelaxation.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/View.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   class ExactHessian;

   NoRelaxation::NoRelaxation(const Model& model, Options& options):
         ConstraintRelaxationStrategy(options),
         original_problem(model),
         globalization_strategy(options) {
   }

   NoRelaxation::~NoRelaxation() = default;

   void NoRelaxation::initialize(Statistics& statistics, Iterate& initial_iterate, bool uses_trust_region,
         EvaluationCache& evaluation_cache, Options& options) {
      this->initial_point.resize(this->original_problem.number_variables);

      // reformulation of the original problem
      this->inequality_handling_method = InequalityHandlingMethodFactory::create(this->original_problem, uses_trust_region,
         1., options);
      //initial_iterate.set_number_variables(this->reformulated_problem->number_variables);

      // initial iterate
      this->inequality_handling_method->evaluate_progress_measures(initial_iterate, evaluation_cache.current_evaluations);
      this->compute_residuals(this->original_problem, initial_iterate, evaluation_cache.current_evaluations);
      this->globalization_strategy.initialize(statistics, initial_iterate);

      // statistics
      this->inequality_handling_method->initialize_statistics(statistics);
   }

   Direction& NoRelaxation::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
         double trust_region_radius, Evaluations& current_evaluations, WarmstartInformation& warmstart_information) {
      DEBUG << "Solving the subproblem\n";
      const bool parameterization_updated = this->inequality_handling_method->update_parameterization(statistics,
         current_iterate);
      // if the problem definition changed, reset the globalization strategy and recompute the current auxiliary measure
      if (parameterization_updated) {
         this->globalization_strategy.reset();
         this->inequality_handling_method->evaluate_progress_measures(current_iterate, current_evaluations); // TODO
      }

      this->initial_point.fill(0.);
      Direction& direction = this->inequality_handling_method->solve(statistics, current_iterate, trust_region_radius,
         this->initial_point, current_evaluations, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, this->original_problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
      warmstart_information.no_changes();
      return direction;
   }

   bool NoRelaxation::solving_feasibility_problem() const {
      return false;
   }

   void NoRelaxation::switch_to_feasibility_problem(Statistics& /*statistics*/, Iterate& /*current_iterate*/,
         Direction& /*direction*/, Evaluations& /*current_evaluations*/, WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("Switching to the feasibility problem should not happen");
   }

   bool NoRelaxation::has_second_order_corrections() const {
      return this->inequality_handling_method->has_second_order_corrections();
   }

   const Direction& NoRelaxation::compute_second_order_correction(Iterate& current_iterate, const Vector<double>& constraints_SOC) {
      return this->inequality_handling_method->compute_second_order_correction(current_iterate, constraints_SOC);
   }

   bool NoRelaxation::is_iterate_acceptable(Statistics& statistics, const Model& /*model*/, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, bool uses_trust_region,
         Evaluations& current_evaluations, Evaluations& trial_evaluations, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) {
      const bool accept_iterate = this->inequality_handling_method->is_iterate_acceptable(statistics, this->globalization_strategy,
         current_iterate, trial_iterate, direction, step_length, current_evaluations, trial_evaluations);
      this->compute_residuals(this->original_problem, trial_iterate, trial_evaluations);
      trial_iterate.status = this->check_termination(this->original_problem, trial_iterate, trial_evaluations);
      if (accept_iterate) {
         user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers,
            this->original_problem.get_objective_multiplier(), trial_iterate.progress.infeasibility,
            trial_iterate.residuals.stationarity, trial_iterate.residuals.complementarity);
      }
      if (uses_trust_region || accept_iterate) {
         this->inequality_handling_method->notify_trial_iterate(statistics, current_iterate, trial_iterate, current_evaluations,
            trial_evaluations);
      }
      warmstart_information.no_changes();
      return accept_iterate;
   }

   std::string NoRelaxation::get_name() const {
      return this->globalization_strategy.get_name() + " " + this->inequality_handling_method->get_name();
   }
} // namespace