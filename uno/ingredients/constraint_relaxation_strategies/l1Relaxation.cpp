// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategyFactory.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"

namespace uno {
   l1Relaxation::l1Relaxation(const Options& options) :
         ConstraintRelaxationStrategy(options),
         inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         inertia_correction_strategy(InertiaCorrectionStrategyFactory::create(options)),
         constraint_violation_coefficient(options.get_double(("l1_constraint_violation_coefficient"))) {
   }

   void l1Relaxation::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate,
         Direction& direction, double trust_region_radius, const Options& options) {
      this->relaxed_problem = std::make_unique<const l1RelaxedProblem>(model, this->penalty_parameter,
         this->constraint_violation_coefficient);
      assert(this->relaxed_problem != nullptr);
      this->feasibility_problem = std::make_unique<const l1RelaxedProblem>(model, 0., 1.);
      assert(this->feasibility_problem != nullptr);

      // Hessian model
      this->hessian_model = HessianModelFactory::create(model, options);

      // memory allocation
      this->inequality_handling_method->initialize(*this->relaxed_problem, initial_iterate, *this->hessian_model,
         *this->inertia_correction_strategy, trust_region_radius);
      direction = Direction(this->relaxed_problem->number_variables, this->relaxed_problem->number_constraints);

      // statistics
      this->inertia_correction_strategy->initialize_statistics(statistics, options);
      this->inequality_handling_method->initialize_statistics(statistics, options);

      // initial iterate
      this->inequality_handling_method->generate_initial_iterate(initial_iterate);
      initial_iterate.evaluate_objective_gradient(model);
      initial_iterate.evaluate_constraints(model);
      this->inequality_handling_method->evaluate_constraint_jacobian(initial_iterate);
      this->relaxed_problem->evaluate_lagrangian_gradient(initial_iterate.residuals.lagrangian_gradient, *this->inequality_handling_method,
         initial_iterate);
      this->compute_primal_dual_residuals(*this->relaxed_problem, initial_iterate);
   }

   void l1Relaxation::compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& /*globalization_strategy*/,
         Iterate& current_iterate, Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information) {
      direction.reset();
      DEBUG << "Solving the subproblem\n";
      direction.set_dimensions(this->relaxed_problem->number_variables, this->relaxed_problem->number_constraints);
      this->inequality_handling_method->solve(statistics, current_iterate, direction, *this->hessian_model,
         *this->inertia_correction_strategy, trust_region_radius, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, this->relaxed_problem->get_number_original_variables()));
      DEBUG3 << direction << '\n';
      warmstart_information.no_changes();
   }

   bool l1Relaxation::solving_feasibility_problem() const {
      return false;
   }

   void l1Relaxation::switch_to_feasibility_problem(Statistics& /*statistics*/, GlobalizationStrategy& /*globalization_strategy*/,
         Iterate& /*current_iterate*/, double /*trust_region_radius*/,
         WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("The problem is unconstrained, switching to the feasibility problem should not happen");
   }

   bool l1Relaxation::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         double trust_region_radius, const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      const bool accept_iterate = this->inequality_handling_method->is_iterate_acceptable(statistics, globalization_strategy,
         *this->hessian_model, *this->inertia_correction_strategy, trust_region_radius, current_iterate, trial_iterate,
         direction, step_length, user_callbacks);
      trial_iterate.status = this->check_termination(model, trial_iterate);
      warmstart_information.no_changes();
      return accept_iterate;
   }

   SolutionStatus l1Relaxation::check_termination(const Model& model, Iterate& iterate) {
      iterate.evaluate_objective_gradient(model);
      iterate.evaluate_constraints(model);

      this->relaxed_problem->evaluate_lagrangian_gradient(iterate.residuals.lagrangian_gradient, *this->inequality_handling_method,
         iterate);
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(*this->relaxed_problem, iterate);
      return ConstraintRelaxationStrategy::check_termination(*this->relaxed_problem, iterate);
   }

   std::string l1Relaxation::get_name() const {
      return this->inequality_handling_method->get_name() + " with " + this->hessian_model->name + " Hessian and " +
         this->inertia_correction_strategy->get_name() + " regularization";
   }

   size_t l1Relaxation::get_number_subproblems_solved() const {
      return this->inequality_handling_method->number_subproblems_solved;
   }
} // namespace