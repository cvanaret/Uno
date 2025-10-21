
#include "InequalityHandlingMethod.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   // infeasibility measure: constraint violation
   void InequalityHandlingMethod::set_infeasibility_measure(const Model& model, Iterate& iterate) const {
      iterate.evaluate_constraints(model);
      iterate.progress.infeasibility = model.constraint_violation(iterate.evaluations.constraints, Norm::L1 /*this->progress_norm*/);
   }

   // objective measure: scaled objective
   void InequalityHandlingMethod::set_objective_measure(const Model& model, Iterate& iterate) const {
      iterate.evaluate_objective(model);
      const double objective = iterate.evaluations.objective;
      iterate.progress.objective = [=](double objective_multiplier) {
         return objective_multiplier * objective;
      };
   }

   void InequalityHandlingMethod::compute_progress_measures(const OptimizationProblem& problem, Iterate& iterate) {
      this->set_infeasibility_measure(problem.model, iterate);
      this->set_objective_measure(problem.model, iterate);
      this->set_auxiliary_measure(iterate);
   }

   bool InequalityHandlingMethod::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const OptimizationProblem& problem, HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy,
         double trust_region_radius, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, UserCallbacks& user_callbacks) {
      this->postprocess_iterate(trial_iterate);
      const double objective_multiplier = problem.get_objective_multiplier();
      trial_iterate.objective_multiplier = objective_multiplier;
      // compute the progress measures
      if (this->subproblem_definition_changed) {
         DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
         globalization_strategy.reset();
         this->set_auxiliary_measure(current_iterate);
         this->subproblem_definition_changed = false;
      }
      this->compute_progress_measures(problem, trial_iterate);

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_iterate.evaluate_objective(problem.model);
         accept_iterate = true;
         statistics.set("status", "0 primal step");
      }
      else {
         const ProgressMeasures predicted_reductions = this->compute_predicted_reductions(problem, hessian_model,
            inertia_correction_strategy, trust_region_radius, current_iterate, direction, step_length);
         accept_iterate = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reductions, objective_multiplier);
      }
      if (accept_iterate) {
         user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers, objective_multiplier,
            trial_iterate.progress.infeasibility, trial_iterate.residuals.stationarity, trial_iterate.residuals.complementarity);
      }
      return accept_iterate;
   }
} // namespace