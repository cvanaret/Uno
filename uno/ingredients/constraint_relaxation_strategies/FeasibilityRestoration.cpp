// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/constraint_relaxation_strategies/l1RelaxedProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategyFactory.hpp"
#include "ingredients/inertia_correction_strategies/UnstableRegularization.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/Expression.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   FeasibilityRestoration::FeasibilityRestoration(const Options& options) :
         ConstraintRelaxationStrategy(options),
         constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")),
         optimality_hessian_model(HessianModelFactory::create(options)),
         feasibility_hessian_model(HessianModelFactory::create(options)),
         optimality_inertia_correction_strategy(InertiaCorrectionStrategyFactory::create(options)),
         feasibility_inertia_correction_strategy(InertiaCorrectionStrategyFactory::create(options)),
         optimality_inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         feasibility_inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         linear_feasibility_tolerance(options.get_double("primal_tolerance")),
         switch_to_optimality_requires_linearized_feasibility(options.get_bool("switch_to_optimality_requires_linearized_feasibility")) {
   }

   void FeasibilityRestoration::initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate,
         Direction& direction, double trust_region_radius, const Options& options) {
      this->optimality_problem = std::make_unique<const OptimizationProblem>(model);
      this->reference_optimality_primals.resize(this->optimality_problem->number_variables);
      this->feasibility_problem = std::make_unique<l1RelaxedProblem>(model, 0., this->constraint_violation_coefficient);

      // memory allocation
      this->optimality_hessian_model->initialize(model);
      this->optimality_inequality_handling_method->initialize(*this->optimality_problem, initial_iterate,
         *this->optimality_hessian_model, *this->optimality_inertia_correction_strategy, trust_region_radius);
      direction = Direction(
         std::max(this->optimality_problem->number_variables, this->feasibility_problem->number_variables),
         std::max(this->optimality_problem->number_constraints, this->feasibility_problem->number_constraints)
      );

      // statistics
      this->optimality_inertia_correction_strategy->initialize_statistics(statistics, options);
      this->optimality_inequality_handling_method->initialize_statistics(statistics, options);
      this->feasibility_inertia_correction_strategy->initialize_statistics(statistics, options);
      this->feasibility_inequality_handling_method->initialize_statistics(statistics, options);
      statistics.add_column("phase", Statistics::int_width, options.get_int("statistics_restoration_phase_column_order"));
      statistics.set("phase", "OPT");

      // initial iterate
      this->optimality_inequality_handling_method->generate_initial_iterate(initial_iterate);
      initial_iterate.evaluate_objective_gradient(model);
      initial_iterate.evaluate_constraints(model);
      this->optimality_inequality_handling_method->evaluate_constraint_jacobian(initial_iterate);
      this->optimality_problem->evaluate_lagrangian_gradient(initial_iterate.residuals.lagrangian_gradient,
         *this->optimality_inequality_handling_method, initial_iterate);
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(*this->optimality_problem, initial_iterate);
   }

   void FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, Direction& direction, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      direction.reset();
      // if we are in the optimality phase, solve the optimality problem
      if (this->current_phase == Phase::OPTIMALITY) {
         statistics.set("phase", "OPT");
         try {
            DEBUG << "Solving the optimality subproblem\n";
            this->solve_subproblem(statistics, *this->optimality_inequality_handling_method, *this->optimality_problem,
               current_iterate, direction, *this->optimality_hessian_model, *this->optimality_inertia_correction_strategy,
               trust_region_radius, warmstart_information);
            if (direction.status == SubproblemStatus::INFEASIBLE) {
               // switch to the feasibility problem, starting from the current direction
               statistics.set("status", std::string("infeasible subproblem"));
               DEBUG << "/!\\ The subproblem is infeasible\n";
               this->switch_to_feasibility_problem(statistics, globalization_strategy, model, current_iterate,
                  trust_region_radius, warmstart_information);
               this->feasibility_inequality_handling_method->set_initial_point(direction.primals);
            }
            else {
               warmstart_information.no_changes();
               return;
            }
         }
         catch (const UnstableRegularization&) {
            this->switch_to_feasibility_problem(statistics, globalization_strategy, model, current_iterate,
               trust_region_radius, warmstart_information);
         }
      }

      // solve the feasibility problem (minimize the constraint violation)
      DEBUG << "Solving the feasibility subproblem\n";
      statistics.set("phase", "FEAS");
      // note: failure of regularization should not happen here, since the feasibility Jacobian has full rank
      this->feasibility_problem->set_proximal_coefficient(this->optimality_inequality_handling_method->proximal_coefficient());
      this->solve_subproblem(statistics, *this->feasibility_inequality_handling_method, *this->feasibility_problem,
         current_iterate, direction, *this->feasibility_hessian_model, *this->feasibility_inertia_correction_strategy,
         trust_region_radius, warmstart_information);
   }

   bool FeasibilityRestoration::solving_feasibility_problem() const {
      return (this->current_phase == Phase::FEASIBILITY_RESTORATION);
   }

   // precondition: this->current_phase == Phase::OPTIMALITY
   void FeasibilityRestoration::switch_to_feasibility_problem(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, double trust_region_radius, WarmstartInformation& warmstart_information) {
      DEBUG << "\nSwitching from optimality to restoration phase\n";
      this->current_phase = Phase::FEASIBILITY_RESTORATION;
      globalization_strategy.notify_switch_to_feasibility(current_iterate.progress);

      // save the current point (progress and primals) upon switching
      this->reference_optimality_progress = current_iterate.progress;
      this->reference_optimality_primals = current_iterate.primals;
      this->feasibility_problem->set_proximal_coefficient(this->optimality_inequality_handling_method->proximal_coefficient());
      this->feasibility_problem->set_proximal_center(this->reference_optimality_primals.data());

      current_iterate.set_number_variables(this->feasibility_problem->number_variables);
      // swap the iterate's multipliers and the feasibility multipliers maintained by the class
      this->other_phase_multipliers.constraints.resize(this->feasibility_problem->number_constraints);
      this->other_phase_multipliers.lower_bounds.resize(this->feasibility_problem->number_variables);
      this->other_phase_multipliers.upper_bounds.resize(this->feasibility_problem->number_variables);
      std::swap(current_iterate.multipliers, this->other_phase_multipliers);

      this->feasibility_inequality_handling_method->initialize_feasibility_problem(current_iterate);
      this->feasibility_inequality_handling_method->set_elastic_variable_values(*this->feasibility_problem, current_iterate);

      DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

      // initialize the feasibility ingredients upon the first switch to feasibility restoration
      if (this->first_switch_to_feasibility) {
         this->feasibility_inequality_handling_method->initialize(*this->feasibility_problem, current_iterate,
            *this->feasibility_hessian_model, *this->feasibility_inertia_correction_strategy, trust_region_radius);
         this->feasibility_hessian_model->initialize(model);

         this->first_switch_to_feasibility = false;
      }

      this->feasibility_inequality_handling_method->evaluate_constraint_jacobian(current_iterate);

      if (Logger::level == INFO) statistics.print_current_line();
      warmstart_information.whole_problem_changed();
   }

   void FeasibilityRestoration::solve_subproblem(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, Iterate& current_iterate, Direction& direction, HessianModel& hessian_model,
         InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      direction.set_dimensions(problem.number_variables, problem.number_constraints);
      inequality_handling_method.solve(statistics, current_iterate, direction, hessian_model, inertia_correction_strategy,
         trust_region_radius, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, problem.get_number_original_variables()));
      DEBUG3 << direction << '\n';
   }

   bool FeasibilityRestoration::can_switch_to_optimality_phase(const Iterate& current_iterate, const GlobalizationStrategy& globalization_strategy,
         const Model& model, const Iterate& trial_iterate, const Direction& direction, double step_length) const {
      if (globalization_strategy.is_infeasibility_sufficiently_reduced(this->reference_optimality_progress, trial_iterate.progress)) {
         if (!this->switch_to_optimality_requires_linearized_feasibility) {
            return true;
         }
         // compute the linearized constraint violation
         // TODO preallocate
         Vector<double> result(model.number_constraints);
         this->feasibility_inequality_handling_method->compute_constraint_jacobian_vector_product(direction.primals, result);
         const double trial_linearized_constraint_violation = model.constraint_violation(current_iterate.evaluations.constraints +
            step_length * result, this->residual_norm);
         return (trial_linearized_constraint_violation <= this->linear_feasibility_tolerance);
      }
      return false;
   }

   void FeasibilityRestoration::switch_back_to_optimality_phase(Iterate& current_iterate, GlobalizationStrategy& globalization_strategy,
         Iterate& trial_iterate) {
      DEBUG << "Switching from restoration back to optimality phase\n";
      this->current_phase = Phase::OPTIMALITY;
      globalization_strategy.notify_switch_to_optimality(current_iterate.progress);

      current_iterate.set_number_variables(this->optimality_problem->number_variables);
      // swap the iterate's multipliers and the optimality multipliers maintained by the class
      std::swap(current_iterate.multipliers, this->other_phase_multipliers);
      trial_iterate.set_number_variables(this->optimality_problem->number_variables);
      current_iterate.objective_multiplier = trial_iterate.objective_multiplier = 1.;
   }

   bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, double trust_region_radius, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      bool accept_iterate = false;
      // determine acceptability, depending on the current phase
      if (this->current_phase == Phase::OPTIMALITY) {
         accept_iterate = this->optimality_inequality_handling_method->is_iterate_acceptable(statistics, globalization_strategy,
            *this->optimality_hessian_model, *this->optimality_inertia_correction_strategy,
            trust_region_radius, current_iterate, trial_iterate, direction, step_length, user_callbacks);
      }
      else {
         accept_iterate = this->feasibility_inequality_handling_method->is_iterate_acceptable(statistics, globalization_strategy,
            *this->feasibility_hessian_model, *this->feasibility_inertia_correction_strategy,
            trust_region_radius, current_iterate, trial_iterate, direction, step_length, user_callbacks);
      }
      trial_iterate.status = this->check_termination(model, trial_iterate);

      // possibly go from restoration phase to optimality phase
      if (trial_iterate.status == SolutionStatus::NOT_OPTIMAL && this->current_phase == Phase::FEASIBILITY_RESTORATION &&
            this->can_switch_to_optimality_phase(current_iterate, globalization_strategy, model, trial_iterate, direction, step_length)) {
         this->switch_back_to_optimality_phase(current_iterate, globalization_strategy, trial_iterate);
         // set a cold start in the subproblem solver
         warmstart_information.whole_problem_changed();
      }
      else {
         warmstart_information.no_changes();
      }
      return accept_iterate;
   }

   SolutionStatus FeasibilityRestoration::check_termination(const Model& model, Iterate& iterate) {
      iterate.evaluate_objective_gradient(model);
      iterate.evaluate_constraints(model);

      if (this->current_phase == Phase::OPTIMALITY) {
         this->optimality_inequality_handling_method->evaluate_constraint_jacobian(iterate);
         this->optimality_problem->evaluate_lagrangian_gradient(iterate.residuals.lagrangian_gradient,
            *this->optimality_inequality_handling_method, iterate);
         ConstraintRelaxationStrategy::compute_primal_dual_residuals(*this->optimality_problem, iterate);
         return ConstraintRelaxationStrategy::check_termination(*this->optimality_problem, iterate);
      }
      else {
         this->feasibility_inequality_handling_method->evaluate_constraint_jacobian(iterate);
         this->feasibility_problem->evaluate_lagrangian_gradient(iterate.residuals.lagrangian_gradient,
            *this->feasibility_inequality_handling_method, iterate);
         ConstraintRelaxationStrategy::compute_primal_dual_residuals(*this->feasibility_problem, iterate);
         return ConstraintRelaxationStrategy::check_termination(*this->feasibility_problem, iterate);
      }
   }

   std::string FeasibilityRestoration::get_name() const {
      return "restoration " + this->optimality_inequality_handling_method->get_name() + " with " +
         this->optimality_hessian_model->get_name() + " Hessian and " + this->optimality_inertia_correction_strategy->get_name() +
         " regularization";
   }

   size_t FeasibilityRestoration::get_hessian_evaluation_count() const {
      return this->optimality_hessian_model->evaluation_count + this->feasibility_hessian_model->evaluation_count;
   }

   size_t FeasibilityRestoration::get_number_subproblems_solved() const {
      return this->optimality_inequality_handling_method->number_subproblems_solved +
         this->feasibility_inequality_handling_method->number_subproblems_solved;
   }
} // namespace