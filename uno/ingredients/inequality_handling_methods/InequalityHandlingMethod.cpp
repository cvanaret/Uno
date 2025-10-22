
#include "InequalityHandlingMethod.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/Options.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Sum.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   InequalityHandlingMethod::InequalityHandlingMethod(const Options& options):
      progress_norm(norm_from_string(options.get_string("progress_norm"))),
      first_order_predicted_reduction(options.get_string("globalization_mechanism") == "LS") {
   }

   void InequalityHandlingMethod::evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate) const {
      problem.set_infeasibility_measure(iterate, this->progress_norm);
      problem.set_objective_measure(iterate);
      problem.set_auxiliary_measure(iterate);
   }

   double InequalityHandlingMethod::compute_predicted_infeasibility_reduction(const Model& model, const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const {
      // predicted infeasibility reduction: "‖c(x)‖ - ‖c(x) + ∇c(x)^T (αd)‖"
      const double current_constraint_violation = model.constraint_violation(current_iterate.evaluations.constraints,
         this->progress_norm);
      Vector<double> result(model.number_constraints);
      this->compute_constraint_jacobian_vector_product(primal_direction, result);
      const double trial_linearized_constraint_violation = model.constraint_violation(current_iterate.evaluations.constraints +
         step_length * result, this->progress_norm);
      return current_constraint_violation - trial_linearized_constraint_violation;
   }

   std::function<double(double)> InequalityHandlingMethod::compute_predicted_objective_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const {
      // predicted objective reduction: "-∇f(x)^T (αd) - α^2/2 d^T H d"
      const double directional_derivative = dot(primal_direction, current_iterate.evaluations.objective_gradient);
      const double quadratic_term = this->first_order_predicted_reduction ? 0. :
         this->compute_hessian_quadratic_product(primal_direction);
      return [=](double objective_multiplier) {
         return step_length * (-objective_multiplier*directional_derivative) - step_length*step_length/2. * quadratic_term;
      };
   }

   ProgressMeasures InequalityHandlingMethod::compute_predicted_reductions(const OptimizationProblem& problem,
         const Iterate& current_iterate, const Direction& direction, double step_length) const {
      return {
         this->compute_predicted_infeasibility_reduction(problem.model, current_iterate, direction.primals, step_length),
         this->compute_predicted_objective_reduction(current_iterate, direction.primals, step_length),
         this->compute_predicted_auxiliary_reduction_model(current_iterate, direction.primals, step_length)
      };
   }

   bool InequalityHandlingMethod::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const OptimizationProblem& problem, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, UserCallbacks& user_callbacks) {
      this->postprocess_iterate(trial_iterate);
      const double objective_multiplier = problem.get_objective_multiplier();

      // evaluate progress measures
      trial_iterate.objective_multiplier = objective_multiplier;
      if (this->subproblem_definition_changed) {
         DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
         globalization_strategy.reset();
         problem.set_auxiliary_measure(current_iterate);
         this->subproblem_definition_changed = false;
      }
      this->evaluate_progress_measures(problem, trial_iterate);

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_iterate.evaluate_objective(problem.model);
         accept_iterate = true;
         statistics.set("status", "0 primal step");
      }
      else {
         const ProgressMeasures predicted_reductions = this->compute_predicted_reductions(problem, current_iterate,
            direction, step_length);
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