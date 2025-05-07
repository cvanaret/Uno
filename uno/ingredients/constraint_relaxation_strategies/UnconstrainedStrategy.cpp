// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "UnconstrainedStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   UnconstrainedStrategy::UnconstrainedStrategy(size_t number_constraints, size_t number_bounds_constraints, const Options& options) :
         ConstraintRelaxationStrategy(number_constraints, number_bounds_constraints, options),
         convexify(options.get_string("inequality_handling_method") != "primal_dual_interior_point" &&
            (options.get_string("globalization_mechanism") != "TR" || options.get_bool("convexify_QP"))),
         hessian_model(HessianModelFactory::create(options.get_string("hessian_model"), this->convexify, options)) {
   }

   void UnconstrainedStrategy::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, Direction& direction, const Options& options) {
      // memory allocation
      this->hessian_model->initialize(model);
      const OptimizationProblem problem{model};
      this->inequality_handling_method->initialize(problem, *this->hessian_model);
      direction = Direction(problem.number_variables, problem.number_constraints);

      // statistics
      this->inequality_handling_method->initialize_statistics(statistics, options);

      // initial iterate
      this->inequality_handling_method->generate_initial_iterate(problem, initial_iterate);
      this->evaluate_progress_measures(model, initial_iterate);
      this->compute_primal_dual_residuals(model, initial_iterate);
      this->set_statistics(statistics, model, initial_iterate);
      this->globalization_strategy->initialize(statistics, initial_iterate, options);
   }

   void UnconstrainedStrategy::compute_feasible_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Direction& direction, WarmstartInformation& warmstart_information) {
      direction.reset();
      DEBUG << "Solving the subproblem\n";
      const OptimizationProblem problem{model};
      this->solve_subproblem(statistics, problem, current_iterate, current_iterate.multipliers, direction,
         *this->hessian_model, warmstart_information);
      warmstart_information.no_changes();
   }

   bool UnconstrainedStrategy::solving_feasibility_problem() const {
      return false;
   }

   void UnconstrainedStrategy::switch_to_feasibility_problem(Statistics& /*statistics*/, const Model& /*model*/, Iterate& /*current_iterate*/,
         WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("The problem is unconstrained, switching to the feasibility problem should not happen");
   }

   void UnconstrainedStrategy::solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, HessianModel& hessian_model, WarmstartInformation& warmstart_information) {
      direction.set_dimensions(problem.number_variables, problem.number_constraints);
      this->inequality_handling_method->solve(statistics, problem, current_iterate, current_multipliers, direction, hessian_model, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
   }

   bool UnconstrainedStrategy::is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      const OptimizationProblem problem{model};
      this->inequality_handling_method->postprocess_iterate(problem, trial_iterate.primals, trial_iterate.multipliers);
      this->compute_progress_measures(model, current_iterate, trial_iterate);
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
         const ProgressMeasures predicted_reduction = this->compute_predicted_reduction_models(model, current_iterate, direction, step_length);
         accept_iterate = this->globalization_strategy->is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
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

   void UnconstrainedStrategy::evaluate_progress_measures(const Model& model, Iterate& iterate) const {
      this->set_infeasibility_measure(model, iterate);
      this->set_objective_measure(model, iterate);
      this->inequality_handling_method->set_auxiliary_measure(model, iterate);
   }

   ProgressMeasures UnconstrainedStrategy::compute_predicted_reduction_models(const Model& model, const Iterate& current_iterate,
         const Direction& direction, double step_length) const {
      return {
         this->compute_predicted_infeasibility_reduction_model(model, current_iterate, direction.primals, step_length),
         this->compute_predicted_objective_reduction_model(current_iterate, direction.primals, step_length),
         this->inequality_handling_method->compute_predicted_auxiliary_reduction_model(model, current_iterate, direction.primals, step_length)
      };
   }

   void UnconstrainedStrategy::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const {
      statistics.set("stationarity", iterate.residuals.stationarity);
      statistics.set("complementarity", iterate.residuals.complementarity);
   }

   std::string UnconstrainedStrategy::get_strategy_combination() const {
      return "unconstrained " + this->globalization_strategy->get_strategy_combination() + " " +
         this->inequality_handling_method->get_strategy_combination();
   }

   size_t UnconstrainedStrategy::get_hessian_evaluation_count() const {
      return this->hessian_model->evaluation_count;
   }
} // namespace