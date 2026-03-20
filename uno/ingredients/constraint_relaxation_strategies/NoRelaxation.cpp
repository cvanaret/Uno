// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <memory>
#include "NoRelaxation.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategyFactory.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/SubproblemSolverFactory.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"

namespace uno {
   NoRelaxation::NoRelaxation(const Model& model, Options& options):
         ConstraintRelaxationStrategy(options),
         original_problem(model),
         inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         hessian_model(HessianModelFactory::create(model, 1., options)),
         inertia_correction_strategy(InertiaCorrectionStrategyFactory::create(options)),
         globalization_strategy(options) {
   }

   void NoRelaxation::initialize(Statistics& statistics, Iterate& initial_iterate, Direction& direction,
         bool uses_trust_region, EvaluationCache& evaluation_cache, const Options& options) {
      // memory allocation
      this->inequality_handling_method->check_problem(this->original_problem, uses_trust_region);
      direction = Direction(this->original_problem.number_variables, this->original_problem.number_constraints);

      // statistics
      this->inertia_correction_strategy->initialize_statistics(statistics);
      this->inequality_handling_method->initialize_statistics(statistics);
      this->hessian_model->initialize_statistics(statistics);

      // initial iterate
      this->reformulated_problem->generate_initial_iterate(initial_iterate, evaluation_cache.current_evaluations);
      this->compute_residuals(this->original_problem, initial_iterate, evaluation_cache.current_evaluations);
      this->globalization_strategy.initialize(statistics, initial_iterate);

      // reformulation
      this->reformulated_problem = this->inequality_handling_method->reformulate(this->original_problem, this->parameterization);
      this->subproblem = std::make_unique<Subproblem>(*this->reformulated_problem, initial_iterate, *this->hessian_model,
         *this->inertia_correction_strategy);
      this->subproblem_solver = SubproblemSolverFactory::create(*this->subproblem, options);
      this->subproblem_solver->initialize_memory(*this->subproblem);
   }

   void NoRelaxation::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double trust_region_radius, Evaluations& current_evaluations, WarmstartInformation& warmstart_information) {
      direction.reset();
      DEBUG << "Solving the subproblem\n";
      direction.set_dimensions(this->original_problem.number_variables, this->original_problem.number_constraints);
      const bool parameterization_updated = this->inequality_handling_method->update_parameterization(statistics,
         this->original_problem, current_iterate, this->parameterization);
      // if the problem definition changed, reset the globalization strategy and recompute the current auxiliary measure
      if (parameterization_updated) {
         this->globalization_strategy.reset();
         this->subproblem->problem.set_auxiliary_measure(current_iterate);
      }
      const Vector<double> initial_point(this->subproblem->number_variables, 0.); // TODO
      this->subproblem_solver->solve(statistics, *this->subproblem, trust_region_radius, initial_point, direction,
         current_evaluations, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, this->original_problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
      warmstart_information.no_changes();
   }

   bool NoRelaxation::solving_feasibility_problem() const {
      return false;
   }

   void NoRelaxation::switch_to_feasibility_problem(Statistics& /*statistics*/, Iterate& /*current_iterate*/,
         Evaluations& /*current_evaluations*/, WarmstartInformation& /*warmstart_information*/) {
      throw std::runtime_error("Switching to the feasibility problem should not happen");
   }

   bool NoRelaxation::is_iterate_acceptable(Statistics& statistics, const Model& /*model*/, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, EvaluationCache& evaluation_cache,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      const bool accept_iterate = ConstraintRelaxationStrategy::is_iterate_acceptable(statistics, this->globalization_strategy,
         *this->subproblem, this->subproblem_solver->get_workspace(), current_iterate, trial_iterate, direction, step_length,
         evaluation_cache, user_callbacks);
      this->compute_residuals(this->original_problem, trial_iterate, evaluation_cache.trial_evaluations);
      trial_iterate.status = this->check_termination(this->original_problem, trial_iterate, evaluation_cache.trial_evaluations);
      if (accept_iterate) {
         this->hessian_model->notify_accepted_iterate(statistics, current_iterate, trial_iterate, evaluation_cache);
      }
      warmstart_information.no_changes();
      return accept_iterate;
   }

   std::string NoRelaxation::get_name() const {
      return this->globalization_strategy.get_name() + " " + this->inequality_handling_method->get_name() + " with " +
         this->hessian_model->name + " Hessian and " + this->inertia_correction_strategy->get_name() + " regularization";
   }

   size_t NoRelaxation::get_number_subproblems_solved() const {
      return this->inequality_handling_method->number_subproblems_solved;
   }
} // namespace