// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "UnconstrainedStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategyFactory.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   UnconstrainedStrategy::UnconstrainedStrategy(size_t number_bound_constraints, const Options& options) :
         ConstraintRelaxationStrategy(options),
         inequality_handling_method(InequalityHandlingMethodFactory::create(number_bound_constraints, options)),
         hessian_model(HessianModelFactory::create(options)),
         regularization_strategy(RegularizationStrategyFactory::create(options)) {
   }

   void UnconstrainedStrategy::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate,
         Direction& direction, const Options& options) {
      const OptimizationProblem problem{model};

      // memory allocation
      this->hessian_model->initialize(model);
      this->inequality_handling_method->initialize(problem, *this->hessian_model, *this->regularization_strategy);
      direction = Direction(problem.number_variables, problem.number_constraints);

      // statistics
      this->regularization_strategy->initialize_statistics(statistics, options);
      this->inequality_handling_method->initialize_statistics(statistics, options);

      // initial iterate
      this->inequality_handling_method->generate_initial_iterate(problem, initial_iterate);
      this->evaluate_progress_measures(*this->inequality_handling_method, model, initial_iterate);
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
         *this->hessian_model, *this->regularization_strategy, trust_region_radius, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
   }

   bool UnconstrainedStrategy::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      const OptimizationProblem problem{model};
      const bool accept_iterate = ConstraintRelaxationStrategy::is_iterate_acceptable(statistics, globalization_strategy,
         model, problem, *this->inequality_handling_method, current_iterate, trial_iterate, trial_iterate.multipliers,
         direction, step_length, user_callbacks);
      ConstraintRelaxationStrategy::set_primal_statistics(statistics, model, trial_iterate);
      warmstart_information.no_changes();
      return accept_iterate;
   }

   void UnconstrainedStrategy::compute_primal_dual_residuals(const Model& model, Iterate& iterate) {
      const OptimizationProblem problem{model};
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(model, problem, problem, iterate);
   }

   void UnconstrainedStrategy::evaluate_progress_measures(InequalityHandlingMethod& inequality_handling_method, const Model& model, Iterate& iterate) const {
      this->set_infeasibility_measure(model, iterate);
      this->set_objective_measure(model, iterate);
      inequality_handling_method.set_auxiliary_measure(model, iterate);
   }

   void UnconstrainedStrategy::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const {
      statistics.set("stationarity", iterate.residuals.stationarity);
      statistics.set("complementarity", iterate.residuals.complementarity);
   }

   std::string UnconstrainedStrategy::get_name() const {
      return this->inequality_handling_method->get_name() + " with " + this->hessian_model->get_name() + " Hessian and " +
         this->regularization_strategy->get_name() + " regularization";
   }

   size_t UnconstrainedStrategy::get_hessian_evaluation_count() const {
      return this->hessian_model->evaluation_count;
   }

   size_t UnconstrainedStrategy::get_number_subproblems_solved() const {
      return this->inequality_handling_method->number_subproblems_solved;
   }
} // namespace