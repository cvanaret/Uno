// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "UnconstrainedStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "layers/SubproblemLayer.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Statistics.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   UnconstrainedStrategy::UnconstrainedStrategy(size_t number_bound_constraints, const Options& options) :
         ConstraintRelaxationStrategy(options),
         inequality_handling_method(InequalityHandlingMethodFactory::create(number_bound_constraints, options)),
         subproblem_layer(options) {
   }

   void UnconstrainedStrategy::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate,
         Direction& direction, const Options& options) {
      const OptimizationProblem problem{model};

      // memory allocation
      this->subproblem_layer.hessian_model->initialize(model);
      this->inequality_handling_method->initialize(problem, *this->subproblem_layer.hessian_model,
         *this->subproblem_layer.regularization_strategy);
      direction = Direction(problem.number_variables, problem.number_constraints);

      // statistics
      this->subproblem_layer.initialize_statistics(statistics, options);
      this->inequality_handling_method->initialize_statistics(statistics, options);

      // initial iterate
      this->inequality_handling_method->generate_initial_iterate(problem, initial_iterate);
      this->evaluate_progress_measures(problem, *this->inequality_handling_method, model, initial_iterate);
      this->compute_primal_dual_residuals(model, initial_iterate);
      this->set_statistics(statistics, model, initial_iterate);
   }

   void UnconstrainedStrategy::compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& /*globalization_strategy*/,
         const Model& model, Iterate& current_iterate, Direction& direction, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      direction.reset();
      DEBUG << "Solving the subproblem\n";
      const OptimizationProblem problem{model};
      this->solve_subproblem(statistics, problem, current_iterate, current_iterate.multipliers, direction, trust_region_radius,
         warmstart_information);
      warmstart_information.no_changes();
   }

   bool UnconstrainedStrategy::solving_feasibility_problem() const {
      return false;
   }

   void UnconstrainedStrategy::switch_to_feasibility_problem(Statistics& /*statistics*/, GlobalizationStrategy& /*globalization_strategy*/,
         const Model& /*model*/, Iterate& /*current_iterate*/, WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("The problem is unconstrained, switching to the feasibility problem should not happen");
   }

   void UnconstrainedStrategy::solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information) {
      direction.set_dimensions(problem.number_variables, problem.number_constraints);
      this->inequality_handling_method->solve(statistics, problem, current_iterate, current_multipliers, direction,
         this->subproblem_layer, trust_region_radius, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
   }

   bool UnconstrainedStrategy::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      const OptimizationProblem problem{model};
      this->inequality_handling_method->postprocess_iterate(problem, trial_iterate.primals, trial_iterate.multipliers);
      this->compute_progress_measures(problem, *this->inequality_handling_method, model, globalization_strategy, current_iterate, trial_iterate);
      constexpr double objective_multiplier = 1.;
      trial_iterate.objective_multiplier = objective_multiplier;
      warmstart_information.no_changes();

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_iterate.evaluate_objective(model);
         accept_iterate = true;
         statistics.set("status", "0 primal step");
      }
      else {
         const ProgressMeasures predicted_reduction = this->compute_predicted_reduction_models(*this->inequality_handling_method,
            model, current_iterate, direction, step_length);
         accept_iterate = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
               predicted_reduction, objective_multiplier);
      }
      ConstraintRelaxationStrategy::set_progress_statistics(statistics, model, trial_iterate);
      if (accept_iterate) {
         user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers, objective_multiplier);
      }
      return accept_iterate;
   }

   void UnconstrainedStrategy::compute_primal_dual_residuals(const Model& model, Iterate& iterate) {
      const OptimizationProblem problem{model};
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(model, problem, problem, iterate);
   }

   void UnconstrainedStrategy::evaluate_progress_measures(const OptimizationProblem& problem, InequalityHandlingMethod& inequality_handling_method,
         const Model& model, Iterate& iterate) const {
      this->set_infeasibility_measure(model, iterate);
      this->set_objective_measure(model, iterate);
      iterate.progress.auxiliary = inequality_handling_method.compute_auxiliary_measure(problem, iterate);
   }

   ProgressMeasures UnconstrainedStrategy::compute_predicted_reduction_models(InequalityHandlingMethod& inequality_handling_method,
         const Model& model, const Iterate& current_iterate, const Direction& direction, double step_length) const {
      return {
         this->compute_predicted_infeasibility_reduction(model, current_iterate, direction.primals, step_length),
         this->compute_predicted_objective_reduction(inequality_handling_method, current_iterate, direction.primals, step_length),
         inequality_handling_method.compute_predicted_auxiliary_reduction_model(model, current_iterate, direction.primals, step_length)
      };
   }

   void UnconstrainedStrategy::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const {
      statistics.set("stationarity", iterate.residuals.stationarity);
      statistics.set("complementarity", iterate.residuals.complementarity);
   }

   std::string UnconstrainedStrategy::get_name() const {
      return "unconstrained";
   }

   size_t UnconstrainedStrategy::get_hessian_evaluation_count() const {
      return this->subproblem_layer.get_hessian_evaluation_count();
   }

   size_t UnconstrainedStrategy::get_number_subproblems_solved() const {
      return this->inequality_handling_method->number_subproblems_solved;
   }
} // namespace