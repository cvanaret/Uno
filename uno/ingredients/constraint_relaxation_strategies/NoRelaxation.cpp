// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "NoRelaxation.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategyFactory.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"

namespace uno {
   NoRelaxation::NoRelaxation(const Model& model, const Options& options):
         ConstraintRelaxationStrategy(options),
         problem(model),
         inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         hessian_model(HessianModelFactory::create(model, options)),
         inertia_correction_strategy(InertiaCorrectionStrategyFactory::create(options)),
         globalization_strategy(options) {
   }

   void NoRelaxation::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate,
         Direction& direction, double trust_region_radius) {
      // memory allocation
      this->inequality_handling_method->initialize(this->problem, initial_iterate, *this->hessian_model,
         *this->inertia_correction_strategy, trust_region_radius);
      direction = Direction(this->problem.number_variables, this->problem.number_constraints);

      // statistics
      this->inertia_correction_strategy->initialize_statistics(statistics);
      this->inequality_handling_method->initialize_statistics(statistics);

      // initial iterate
      this->inequality_handling_method->generate_initial_iterate(initial_iterate);
      initial_iterate.evaluate_objective_gradient(model);
      initial_iterate.evaluate_constraints(model);
      this->inequality_handling_method->evaluate_constraint_jacobian(initial_iterate);
      const auto& evaluation_space = this->inequality_handling_method->get_evaluation_space();
      this->problem.evaluate_lagrangian_gradient(initial_iterate.residuals.lagrangian_gradient, evaluation_space, initial_iterate);
      this->compute_primal_dual_residuals(this->problem, initial_iterate);
      this->globalization_strategy.initialize(statistics, initial_iterate);
   }

   void NoRelaxation::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double trust_region_radius, WarmstartInformation& warmstart_information) {
      direction.reset();
      DEBUG << "Solving the subproblem\n";
      direction.set_dimensions(this->problem.number_variables, this->problem.number_constraints);
      this->inequality_handling_method->solve(statistics, current_iterate, direction, trust_region_radius, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, this->problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
      warmstart_information.no_changes();
   }

   bool NoRelaxation::solving_feasibility_problem() const {
      return false;
   }

   void NoRelaxation::switch_to_feasibility_problem(Statistics& /*statistics*/, Iterate& /*current_iterate*/,
         double /*trust_region_radius*/, WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("Switching to the feasibility problem should not happen");
   }

   bool NoRelaxation::is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) {
      const bool accept_iterate = this->inequality_handling_method->is_iterate_acceptable(statistics, this->globalization_strategy,
         current_iterate, trial_iterate, direction, step_length, user_callbacks);
      trial_iterate.status = this->check_termination(model, trial_iterate);
      warmstart_information.no_changes();
      return accept_iterate;
   }

   SolutionStatus NoRelaxation::check_termination(const Model& model, Iterate& iterate) {
      iterate.evaluate_objective_gradient(model);
      iterate.evaluate_constraints(model);

      const auto& evaluation_space = this->inequality_handling_method->get_evaluation_space();
      this->problem.evaluate_lagrangian_gradient(iterate.residuals.lagrangian_gradient, evaluation_space, iterate);
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(this->problem, iterate);
      return ConstraintRelaxationStrategy::check_termination(this->problem, iterate);
   }

   std::string NoRelaxation::get_name() const {
      return this->globalization_strategy.get_name() + " " + this->inequality_handling_method->get_name() + " with " +
         this->hessian_model->name + " Hessian and " + this->inertia_correction_strategy->get_name() + " regularization";
   }

   size_t NoRelaxation::get_number_subproblems_solved() const {
      return this->inequality_handling_method->number_subproblems_solved;
   }
} // namespace