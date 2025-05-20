// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "OptimizationProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Model.hpp"
#include "optimization/Multipliers.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"
#include "symbolic/Expression.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Options& options):
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

   // infeasibility measure: constraint violation
   void ConstraintRelaxationStrategy::set_infeasibility_measure(const Model& model, Iterate& iterate) const {
      iterate.evaluate_constraints(model);
      iterate.progress.infeasibility = model.constraint_violation(iterate.model_evaluations.constraints, this->progress_norm);
   }

   // objective measure: scaled objective
   void ConstraintRelaxationStrategy::set_objective_measure(const Model& model, Iterate& iterate) const {
      iterate.evaluate_objective(model);
      const double objective = iterate.model_evaluations.objective;
      iterate.progress.objective = [=](double objective_multiplier) {
         return objective_multiplier * objective;
      };
   }

   double ConstraintRelaxationStrategy::compute_predicted_infeasibility_reduction(const Model& model, const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const {
      // predicted infeasibility reduction: "‖c(x)‖ - ‖c(x) + ∇c(x)^T (αd)‖"
      const double current_constraint_violation = model.constraint_violation(current_iterate.model_evaluations.constraints, this->progress_norm);
      const double trial_linearized_constraint_violation = model.constraint_violation(current_iterate.model_evaluations.constraints +
         step_length * (current_iterate.model_evaluations.constraint_jacobian * primal_direction), this->progress_norm);
      return current_constraint_violation - trial_linearized_constraint_violation;
   }

   std::function<double(double)> ConstraintRelaxationStrategy::compute_predicted_objective_reduction(InequalityHandlingMethod& inequality_handling_method,
         const Iterate& current_iterate, const Vector<double>& primal_direction, double step_length) const {
      // predicted objective reduction: "-∇f(x)^T (αd) - α^2/2 d^T H d"
      const double directional_derivative = dot(primal_direction, current_iterate.model_evaluations.objective_gradient);
      const double quadratic_term = this->first_order_predicted_reduction ? 0. :
         inequality_handling_method.hessian_quadratic_product(primal_direction);
      return [=](double objective_multiplier) {
         return step_length * (-objective_multiplier*directional_derivative) - step_length*step_length/2. * quadratic_term;
      };
   }

   void ConstraintRelaxationStrategy::compute_progress_measures(const OptimizationProblem& problem, InequalityHandlingMethod& inequality_handling_method,
         const Model& model, GlobalizationStrategy& globalization_strategy, Iterate& current_iterate, Iterate& trial_iterate) {
      if (inequality_handling_method.subproblem_definition_changed) {
         DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
         globalization_strategy.reset();
         current_iterate.progress.auxiliary = inequality_handling_method.compute_auxiliary_measure(problem, current_iterate);
         inequality_handling_method.subproblem_definition_changed = false;
      }
      this->evaluate_progress_measures(problem, inequality_handling_method, model, trial_iterate);
   }

   void ConstraintRelaxationStrategy::compute_primal_dual_residuals(const Model& model, const OptimizationProblem& problem,
         Iterate& iterate) {
      iterate.evaluate_objective_gradient(model);
      iterate.evaluate_constraints(model);
      iterate.evaluate_constraint_jacobian(model);
      
      problem.evaluate_lagrangian_gradient(iterate.residuals.lagrangian_gradient, iterate, iterate.multipliers);
      iterate.residuals.stationarity = OptimizationProblem::stationarity_error(iterate.residuals.lagrangian_gradient, iterate.objective_multiplier,
            this->residual_norm);

      // constraint violation of the original problem
      iterate.primal_feasibility = model.constraint_violation(iterate.model_evaluations.constraints, this->residual_norm);

      // complementarity error
      const double shift_value = 0.;
      iterate.residuals.complementarity = problem.complementarity_error(iterate.primals, iterate.model_evaluations.constraints,
            iterate.multipliers, shift_value, this->residual_norm);

      // scaling factors
      iterate.residuals.stationarity_scaling = this->compute_stationarity_scaling(model, iterate.multipliers);
      iterate.residuals.complementarity_scaling = this->compute_complementarity_scaling(model, iterate.multipliers);
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

   void ConstraintRelaxationStrategy::assemble_trial_iterate(Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double primal_step_length, double dual_step_length) {
      //assert(current_iterate.number_variables == direction.number_variables && "The primal iterate and the primal direction do not match");
      current_iterate.set_number_variables(direction.number_variables);
      trial_iterate.set_number_variables(direction.number_variables);

      // resize the trial iterate to the same size as the current iterate
      trial_iterate.primals.resize(current_iterate.primals.size());
      trial_iterate.multipliers.constraints.resize(current_iterate.multipliers.constraints.size());
      trial_iterate.multipliers.lower_bounds.resize(current_iterate.multipliers.lower_bounds.size());
      trial_iterate.multipliers.upper_bounds.resize(current_iterate.multipliers.upper_bounds.size());
      trial_iterate.residuals.lagrangian_gradient.resize(current_iterate.residuals.lagrangian_gradient.size());

      // take primal step
      trial_iterate.primals = current_iterate.primals + primal_step_length * direction.primals;

      // take dual step: line-search carried out only on constraint multipliers. Bound multipliers updated with full step
      trial_iterate.multipliers.constraints = current_iterate.multipliers.constraints + dual_step_length * direction.multipliers.constraints;
      trial_iterate.multipliers.lower_bounds = current_iterate.multipliers.lower_bounds + direction.multipliers.lower_bounds;
      trial_iterate.multipliers.upper_bounds = current_iterate.multipliers.upper_bounds + direction.multipliers.upper_bounds;

      trial_iterate.progress.reset();
      trial_iterate.is_objective_computed = false;
      trial_iterate.is_objective_gradient_computed = false;
      trial_iterate.are_constraints_computed = false;
      trial_iterate.is_constraint_jacobian_computed = false;
      trial_iterate.status = IterateStatus::NOT_OPTIMAL;
   }

   IterateStatus ConstraintRelaxationStrategy::check_termination(const Model& model, Iterate& iterate) {
      if (iterate.is_objective_computed && iterate.model_evaluations.objective < this->unbounded_objective_threshold) {
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

   void ConstraintRelaxationStrategy::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const {
      statistics.set("stationarity", iterate.residuals.stationarity);
      statistics.set("complementarity", iterate.residuals.complementarity);
   }

   IterateStatus ConstraintRelaxationStrategy::check_first_order_convergence(const Model& model, Iterate& current_iterate, double tolerance) const {
      // evaluate termination conditions based on optimality conditions
      const bool stationarity = (current_iterate.residuals.stationarity / current_iterate.residuals.stationarity_scaling <= tolerance);
      const bool primal_feasibility = (current_iterate.primal_feasibility <= tolerance);
      const bool complementarity = (current_iterate.residuals.complementarity / current_iterate.residuals.complementarity_scaling <= tolerance);
      DEBUG << "\nTermination criteria for tolerance = " << tolerance << ":\n";
      DEBUG << "Stationarity: " << std::boolalpha << stationarity << '\n';
      DEBUG << "Primal feasibility: " << std::boolalpha << primal_feasibility << '\n';
      DEBUG << "Complementarity: " << std::boolalpha << complementarity << '\n';

      if (stationarity && primal_feasibility && 0. < current_iterate.objective_multiplier && complementarity) {
         // feasible regular stationary point
         return IterateStatus::FEASIBLE_KKT_POINT;
      }
      /*
      else if (model.is_constrained() && feasibility_stationarity && !primal_feasibility && feasibility_complementarity && no_trivial_duals) {
         // no primal feasibility, stationary point of constraint violation
         return IterateStatus::INFEASIBLE_STATIONARY_POINT;
      }
      */
      return IterateStatus::NOT_OPTIMAL;
   }

   void ConstraintRelaxationStrategy::set_statistics(Statistics& statistics, const Model& model, const Iterate& iterate) const {
      ConstraintRelaxationStrategy::set_progress_statistics(statistics, model, iterate);
      this->set_dual_residuals_statistics(statistics, iterate);
   }

   void ConstraintRelaxationStrategy::set_progress_statistics(Statistics& statistics, const Model& model, const Iterate& iterate) {
      statistics.set("objective", iterate.model_evaluations.objective);
      if (model.is_constrained()) {
         statistics.set("primal feas", iterate.progress.infeasibility);
      }
   }

   void ConstraintRelaxationStrategy::check_unboundedness(const Direction& direction) {
      if (direction.status == SubproblemStatus::UNBOUNDED_PROBLEM) {
         throw std::runtime_error("The subproblem is unbounded, this should not happen. If the subproblem has curvature,"
            "use regularization. If not, use a trust-region method.\n");
      }
   }
} // namespace