// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "UnconstrainedStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategyFactory.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"

namespace uno {
   UnconstrainedStrategy::UnconstrainedStrategy(const Options& options) :
         ConstraintRelaxationStrategy(options),
         inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         hessian_model(HessianModelFactory::create(options)),
         regularization_strategy(RegularizationStrategyFactory::create(options)) {
   }

   void UnconstrainedStrategy::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate,
         Direction& direction, double trust_region_radius, const Options& options) {
      const OptimizationProblem problem{model};

      // memory allocation
      this->hessian_model->initialize(model);
      this->inequality_handling_method->initialize(problem, initial_iterate, *this->hessian_model,
         *this->regularization_strategy, trust_region_radius);
      direction = Direction(problem.number_variables, problem.number_constraints);

      // statistics
      this->regularization_strategy->initialize_statistics(statistics, options);
      this->inequality_handling_method->initialize_statistics(statistics, options);

      // initial iterate
      this->inequality_handling_method->generate_initial_iterate(problem, initial_iterate);
      this->evaluate_progress_measures(*this->inequality_handling_method, problem, initial_iterate);
      initial_iterate.evaluate_objective_gradient(model);
      initial_iterate.evaluate_constraints(model);
      this->inequality_handling_method->evaluate_constraint_jacobian(problem, initial_iterate);
      problem.evaluate_lagrangian_gradient(initial_iterate.residuals.lagrangian_gradient, *this->inequality_handling_method,
         initial_iterate);
      this->compute_primal_dual_residuals(problem, initial_iterate);
   }

   void UnconstrainedStrategy::compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& /*globalization_strategy*/,
         const Model& model, Iterate& current_iterate, Direction& direction, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      direction.reset();
      DEBUG << "Solving the subproblem\n";
      const OptimizationProblem problem{model};
      direction.set_dimensions(problem.number_variables, problem.number_constraints);
      this->inequality_handling_method->solve(statistics, problem, current_iterate, direction, *this->hessian_model,
         *this->regularization_strategy, trust_region_radius, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
      warmstart_information.no_changes();
   }

   bool UnconstrainedStrategy::solving_feasibility_problem() const {
      return false;
   }

   void UnconstrainedStrategy::switch_to_feasibility_problem(Statistics& /*statistics*/, GlobalizationStrategy& /*globalization_strategy*/,
         const Model& /*model*/, Iterate& /*current_iterate*/, double /*trust_region_radius*/,
         WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("The problem is unconstrained, switching to the feasibility problem should not happen");
   }

   bool UnconstrainedStrategy::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      const OptimizationProblem problem{model};
      const bool accept_iterate = ConstraintRelaxationStrategy::is_iterate_acceptable(statistics, globalization_strategy,
         problem, *this->inequality_handling_method, current_iterate, trial_iterate, direction, step_length, user_callbacks);
      trial_iterate.status = this->check_termination(model, trial_iterate);
      warmstart_information.no_changes();
      return accept_iterate;
   }

   SolutionStatus UnconstrainedStrategy::check_termination(const Model& model, Iterate& iterate) {
      iterate.evaluate_objective_gradient(model);
      iterate.evaluate_constraints(model);

      const OptimizationProblem problem{model};
      problem.evaluate_lagrangian_gradient(iterate.residuals.lagrangian_gradient, *this->inequality_handling_method, iterate);
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(problem, iterate);
      return ConstraintRelaxationStrategy::check_termination(problem, iterate);
   }

   void UnconstrainedStrategy::evaluate_progress_measures(InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, Iterate& iterate) const {
      this->set_infeasibility_measure(problem.model, iterate);
      this->set_objective_measure(problem.model, iterate);
      inequality_handling_method.set_auxiliary_measure(problem, iterate);
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