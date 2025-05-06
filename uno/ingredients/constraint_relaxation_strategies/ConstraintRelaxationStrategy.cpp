// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "OptimizationProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"
#include "optimization/Direction.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Multipliers.hpp"
#include "symbolic/VectorView.hpp"
#include "symbolic/Expression.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Options& options):
         globalization_strategy(GlobalizationStrategyFactory::create(options.get_string("globalization_strategy"), options)),
         inequality_handling_method(InequalityHandlingMethodFactory::create(options)),
         progress_norm(norm_from_string(options.get_string("progress_norm"))),
         residual_norm(norm_from_string(options.get_string("residual_norm"))),
         residual_scaling_threshold(options.get_double("residual_scaling_threshold")),
         tight_tolerance(options.get_double("tolerance")),
         loose_tolerance(options.get_double("loose_tolerance")),
         loose_tolerance_consecutive_iteration_threshold(options.get_unsigned_int("loose_tolerance_consecutive_iteration_threshold")),
         unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")),
         first_order_predicted_reduction(options.get_string("globalization_mechanism") == "LS") {
   }

   ConstraintRelaxationStrategy::~ConstraintRelaxationStrategy() { }

   void ConstraintRelaxationStrategy::set_trust_region_radius(double trust_region_radius) {
      this->inequality_handling_method->set_trust_region_radius(trust_region_radius);
   }

   // with initial point
   void ConstraintRelaxationStrategy::compute_feasible_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
      Direction& direction, const Vector<double>& initial_point, WarmstartInformation& warmstart_information) {
      this->inequality_handling_method->set_initial_point(initial_point);
      this->compute_feasible_direction(statistics, model, current_iterate, direction, warmstart_information);
   }

   // infeasibility measure: constraint violation
   void ConstraintRelaxationStrategy::set_infeasibility_measure(const Model& model, Iterate& iterate) const {
      iterate.evaluate_constraints(model);
      iterate.progress.infeasibility = model.constraint_violation(iterate.evaluations.constraints, this->progress_norm);
   }

   // objective measure: scaled objective
   void ConstraintRelaxationStrategy::set_objective_measure(const Model& model, Iterate& iterate) const {
      iterate.evaluate_objective(model);
      const double objective = iterate.evaluations.objective;
      iterate.progress.objective = [=](double objective_multiplier) {
         return objective_multiplier * objective;
      };
   }

   double ConstraintRelaxationStrategy::compute_predicted_infeasibility_reduction_model(const Model& model, const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const {
      // predicted infeasibility reduction: "‖c(x)‖ - ‖c(x) + ∇c(x)^T (αd)‖"
      const double current_constraint_violation = model.constraint_violation(current_iterate.evaluations.constraints, this->progress_norm);
      const double trial_linearized_constraint_violation = model.constraint_violation(current_iterate.evaluations.constraints +
         step_length * (current_iterate.evaluations.constraint_jacobian * primal_direction), this->progress_norm);
      return current_constraint_violation - trial_linearized_constraint_violation;
   }

   std::function<double(double)> ConstraintRelaxationStrategy::compute_predicted_objective_reduction_model(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const {
      // predicted objective reduction: "-∇f(x)^T (αd) - α^2/2 d^T H d"
      const double directional_derivative = dot(primal_direction, current_iterate.evaluations.objective_gradient);
      const double quadratic_term = this->first_order_predicted_reduction ? 0. :
         this->inequality_handling_method->hessian_quadratic_product(primal_direction);
      return [=](double objective_multiplier) {
         return step_length * (-objective_multiplier*directional_derivative) - step_length*step_length/2. * quadratic_term;
      };
   }

   void ConstraintRelaxationStrategy::compute_progress_measures(const Model& model, Iterate& current_iterate, Iterate& trial_iterate) {
      if (this->inequality_handling_method->subproblem_definition_changed) {
         DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
         this->globalization_strategy->reset();
         this->inequality_handling_method->set_auxiliary_measure(model, current_iterate);
         this->inequality_handling_method->subproblem_definition_changed = false;
      }
      this->evaluate_progress_measures(model, trial_iterate);
   }

   void ConstraintRelaxationStrategy::compute_primal_dual_residuals(const Model& model, const OptimizationProblem& optimality_problem,
         const OptimizationProblem& feasibility_problem, Iterate& iterate) {
      iterate.evaluate_objective_gradient(model);
      iterate.evaluate_constraints(model);
      iterate.evaluate_constraint_jacobian(model);

      // stationarity errors:
      // - for KKT conditions: with standard multipliers and current objective multiplier
      // - for FJ conditions: with standard multipliers and 0 objective multiplier
      // - for feasibility problem: with feasibility multipliers and 0 objective multiplier
      optimality_problem.evaluate_lagrangian_gradient(iterate.residuals.lagrangian_gradient, iterate, iterate.multipliers);
      iterate.residuals.stationarity = OptimizationProblem::stationarity_error(iterate.residuals.lagrangian_gradient, iterate.objective_multiplier,
            this->residual_norm);
      feasibility_problem.evaluate_lagrangian_gradient(iterate.feasibility_residuals.lagrangian_gradient, iterate, iterate.feasibility_multipliers);
      iterate.feasibility_residuals.stationarity = OptimizationProblem::stationarity_error(iterate.feasibility_residuals.lagrangian_gradient, 0.,
            this->residual_norm);

      // constraint violation of the original problem
      iterate.primal_feasibility = model.constraint_violation(iterate.evaluations.constraints, this->residual_norm);

      // complementarity error
      const double shift_value = 0.;
      iterate.residuals.complementarity = optimality_problem.complementarity_error(iterate.primals, iterate.evaluations.constraints,
            iterate.multipliers, shift_value, this->residual_norm);
      iterate.feasibility_residuals.complementarity = feasibility_problem.complementarity_error(iterate.primals, iterate.evaluations.constraints,
            iterate.feasibility_multipliers, shift_value, this->residual_norm);

      // scaling factors
      iterate.residuals.stationarity_scaling = this->compute_stationarity_scaling(model, iterate.multipliers);
      iterate.residuals.complementarity_scaling = this->compute_complementarity_scaling(model, iterate.multipliers);
      iterate.feasibility_residuals.stationarity_scaling = this->compute_stationarity_scaling(model, iterate.feasibility_multipliers);
      iterate.feasibility_residuals.complementarity_scaling = this->compute_complementarity_scaling(model, iterate.feasibility_multipliers);
   }

   double ConstraintRelaxationStrategy::compute_stationarity_scaling(const Model& model, const Multipliers& multipliers) const {
      const size_t total_size = model.get_lower_bounded_variables().size() + model.get_upper_bounded_variables().size() + model.number_constraints;
      if (total_size == 0) {
         return 1.;
      }
      else {
         const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
         const double multiplier_norm = norm_1(
               view(multipliers.constraints, 0, model.number_constraints),
               view(multipliers.lower_bounds, 0, model.number_variables),
               view(multipliers.upper_bounds, 0, model.number_variables)
         );
         return std::max(1., multiplier_norm / scaling_factor);
      }
   }

   double ConstraintRelaxationStrategy::compute_complementarity_scaling(const Model& model, const Multipliers& multipliers) const {
      const size_t total_size = model.get_lower_bounded_variables().size() + model.get_upper_bounded_variables().size();
      if (total_size == 0) {
         return 1.;
      }
      else {
         const double scaling_factor = this->residual_scaling_threshold * static_cast<double>(total_size);
         const double bound_multiplier_norm = norm_1(
               view(multipliers.lower_bounds, 0, model.number_variables),
               view(multipliers.upper_bounds, 0, model.number_variables)
         );
         return std::max(1., bound_multiplier_norm / scaling_factor);
      }
   }

   IterateStatus ConstraintRelaxationStrategy::check_termination(const Model& model, Iterate& iterate) {
      if (iterate.is_objective_computed && iterate.evaluations.objective < this->unbounded_objective_threshold) {
         return IterateStatus::UNBOUNDED;
      }

      // compute the residuals
      this->compute_primal_dual_residuals(model, iterate);

      // test convergence wrt the tight tolerance
      const IterateStatus status_tight_tolerance = this->check_first_order_convergence(model, iterate, this->tight_tolerance);
      if (status_tight_tolerance != IterateStatus::NOT_OPTIMAL || this->loose_tolerance <= this->tight_tolerance) {
         return status_tight_tolerance;
      }

      // if not converged, check convergence wrt loose tolerance (provided it is strictly looser than the tight tolerance)
      const IterateStatus status_loose_tolerance = this->check_first_order_convergence(model, iterate, this->loose_tolerance);
      // if converged, keep track of the number of consecutive iterations
      if (status_loose_tolerance != IterateStatus::NOT_OPTIMAL) {
         this->loose_tolerance_consecutive_iterations++;
      }
      else {
         this->loose_tolerance_consecutive_iterations = 0;
         return IterateStatus::NOT_OPTIMAL;
      }
      // check if loose tolerance achieved for enough consecutive iterations
      if (this->loose_tolerance_consecutive_iteration_threshold <= this->loose_tolerance_consecutive_iterations) {
         return status_loose_tolerance;
      }
      else {
         return IterateStatus::NOT_OPTIMAL;
      }
   }

   IterateStatus ConstraintRelaxationStrategy::check_first_order_convergence(const Model& model, Iterate& current_iterate, double tolerance) const {
      // evaluate termination conditions based on optimality conditions
      const bool stationarity = (current_iterate.residuals.stationarity / current_iterate.residuals.stationarity_scaling <= tolerance);
      const bool primal_feasibility = (current_iterate.primal_feasibility <= tolerance);
      const bool complementarity = (current_iterate.residuals.complementarity / current_iterate.residuals.complementarity_scaling <= tolerance);

      const bool feasibility_stationarity = (current_iterate.feasibility_residuals.stationarity <= tolerance);
      const bool feasibility_complementarity = (current_iterate.feasibility_residuals.complementarity <= tolerance);
      const bool no_trivial_duals = current_iterate.feasibility_multipliers.not_all_zero(model.number_variables, tolerance);

      DEBUG << "\nTermination criteria for tolerance = " << tolerance << ":\n";
      DEBUG << "Stationarity: " << std::boolalpha << stationarity << '\n';
      DEBUG << "Primal feasibility: " << std::boolalpha << primal_feasibility << '\n';
      DEBUG << "Complementarity: " << std::boolalpha << complementarity << '\n';

      DEBUG << "Feasibility stationarity: " << std::boolalpha << feasibility_stationarity << '\n';
      DEBUG << "Feasibility complementarity: " << std::boolalpha << feasibility_complementarity << '\n';
      DEBUG << "Not all zero multipliers: " << std::boolalpha << no_trivial_duals << "\n\n";

      if (stationarity && primal_feasibility && 0. < current_iterate.objective_multiplier && complementarity) {
         // feasible regular stationary point
         return IterateStatus::FEASIBLE_KKT_POINT;
      }
      else if (model.is_constrained() && feasibility_stationarity && !primal_feasibility && feasibility_complementarity && no_trivial_duals) {
         // no primal feasibility, stationary point of constraint violation
         return IterateStatus::INFEASIBLE_STATIONARY_POINT;
      }
      return IterateStatus::NOT_OPTIMAL;
   }

   void ConstraintRelaxationStrategy::set_statistics(Statistics& statistics, const Model& model, const Iterate& iterate) const {
      this->set_progress_statistics(statistics, model, iterate);
      this->set_dual_residuals_statistics(statistics, iterate);
   }

   void ConstraintRelaxationStrategy::set_progress_statistics(Statistics& statistics, const Model& model, const Iterate& iterate) const {
      statistics.set("objective", iterate.evaluations.objective);
      if (model.is_constrained()) {
         statistics.set("primal feas", iterate.progress.infeasibility);
      }
   }

   size_t ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
      return this->inequality_handling_method->number_subproblems_solved;
   }
} // namespace