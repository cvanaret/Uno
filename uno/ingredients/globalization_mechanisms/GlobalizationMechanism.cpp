// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <limits>
#include "GlobalizationMechanism.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Sum.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   GlobalizationMechanism::GlobalizationMechanism(const Model& model, bool use_trust_region, const Options& options,
      std::vector<OptionOverride>& option_overrides):
         constraint_relaxation_strategy(ConstraintRelaxationStrategyFactory::create(model, use_trust_region, options,
            option_overrides)),
         primal_tolerance(options.get_double("primal_tolerance")),
         tight_unscaled_primal_infeasibility_tolerance(options.get_double("tight_unscaled_primal_infeasibility_tolerance")),
         residual_norm(norm_from_string(options.get_string("residual_norm"))) {
   }

   GlobalizationMechanism::~GlobalizationMechanism() = default;

   size_t GlobalizationMechanism::get_number_subproblems_solved() const {
      return this->constraint_relaxation_strategy->get_number_subproblems_solved();
   }

   namespace {
      // the components of vector2 differ from those of vector1 by round-off only
      bool differ_by_roundoff(const Vector<double>& vector1, const Vector<double>& vector2, double factor) {
         if (vector1.size() != vector2.size()) {
            return false;
         }
         constexpr double machine_epsilon = std::numeric_limits<double>::epsilon();
         for (size_t index: Range(vector1.size())) {
            // written so that a NaN is never a round-off difference
            if (!(std::abs(vector2[index] - vector1[index]) <= factor * machine_epsilon * (1. + std::abs(vector1[index])))) {
               return false;
            }
         }
         return true;
      }

      // same primal-dual point up to round-off
      bool differ_by_roundoff(const Iterate& iterate1, const Iterate& iterate2, double factor) {
         return differ_by_roundoff(iterate1.primals, iterate2.primals, factor) &&
            differ_by_roundoff(iterate1.multipliers.constraints, iterate2.multipliers.constraints, factor) &&
            differ_by_roundoff(iterate1.multipliers.lower_bounds, iterate2.multipliers.lower_bounds, factor) &&
            differ_by_roundoff(iterate1.multipliers.upper_bounds, iterate2.multipliers.upper_bounds, factor);
      }
   } // namespace

   void GlobalizationMechanism::detect_stagnation(const Model& model, const Iterate& current_iterate, Iterate& trial_iterate,
         Evaluations& trial_evaluations) {
      // an update of the parameterization (e.g. a decrease of the barrier parameter) is progress
      const size_t number_parameterization_updates = this->constraint_relaxation_strategy->get_number_parameterization_updates();
      const bool parameterization_updated = (number_parameterization_updates != this->stagnation_number_parameterization_updates);
      this->stagnation_number_parameterization_updates = number_parameterization_updates;
      if (trial_iterate.status != SolutionStatus::NOT_OPTIMAL || parameterization_updated ||
            !differ_by_roundoff(current_iterate, trial_iterate, GlobalizationMechanism::stagnation_roundoff_factor)) {
         this->number_stagnant_iterations = 0;
         return;
      }
      ++this->number_stagnant_iterations;
      DEBUG << "The trial iterate differs from the current iterate by round-off only (" << this->number_stagnant_iterations <<
         " consecutive iteration(s))\n";
      if (GlobalizationMechanism::max_stagnant_iterations <= this->number_stagnant_iterations) {
         DEBUG << "Uno stagnates: terminating\n";
         trial_iterate.status = this->small_step_status(model, trial_iterate, trial_evaluations);
      }
   }

   SolutionStatus GlobalizationMechanism::small_step_status(const Model& model, const Iterate& iterate,
         Evaluations& evaluations) const {
      // an iterate that satisfies the acceptable conditions is reported as such (as with acceptable termination)
      const SolutionStatus acceptable_status = this->constraint_relaxation_strategy->get_last_iterate_acceptable_status();
      if (acceptable_status != SolutionStatus::NOT_OPTIMAL) {
         return acceptable_status;
      }
      if (iterate.scaled_residuals.primal_infeasibility <= this->primal_tolerance) {
         return SolutionStatus::FEASIBLE_SMALL_STEP;
      }
      evaluations.evaluate_constraints(model, iterate.primals); // no-op if cached or unconstrained
      const double unscaled_primal_infeasibility = model.unscaled_primal_infeasibility(
         iterate.scaled_residuals.primal_infeasibility, evaluations.constraints, this->residual_norm);
      return (unscaled_primal_infeasibility <= this->tight_unscaled_primal_infeasibility_tolerance) ?
         SolutionStatus::FEASIBLE_SMALL_STEP : SolutionStatus::INFEASIBLE_SMALL_STEP;
   }

   void GlobalizationMechanism::assemble_trial_iterate(Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double primal_step_length, double constraint_dual_step_length, double bound_dual_step_length) {
      trial_iterate.set_dimensions(current_iterate.primals.size(), current_iterate.multipliers.constraints.size());

      // take primal step
      trial_iterate.primals = current_iterate.primals + primal_step_length * direction.primals;
      // take dual step: line-search carried out only on constraint multipliers. Bound multipliers updated with full step
      trial_iterate.multipliers.constraints = current_iterate.multipliers.constraints + constraint_dual_step_length * direction.multipliers.constraints;
      trial_iterate.multipliers.lower_bounds = current_iterate.multipliers.lower_bounds + bound_dual_step_length * direction.multipliers.lower_bounds;
      trial_iterate.multipliers.upper_bounds = current_iterate.multipliers.upper_bounds + bound_dual_step_length * direction.multipliers.upper_bounds;
      trial_iterate.progress.reset();
      trial_iterate.status = SolutionStatus::NOT_OPTIMAL;
   }

   void GlobalizationMechanism::set_primal_statistics(Statistics& statistics, const Model& model, const Iterate& iterate,
         const Evaluations& evaluations) {
      if (evaluations.is_objective_computed) {
         statistics.set("Objective", evaluations.objective);
      }
      if (model.is_constrained()) {
         statistics.set("Infeas", iterate.scaled_residuals.primal_infeasibility);
      }
   }

   void GlobalizationMechanism::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) {
      statistics.set("Statio", iterate.scaled_residuals.stationarity / iterate.scaled_residuals.stationarity_scaling);
      statistics.set("Compl", iterate.scaled_residuals.complementarity / iterate.scaled_residuals.complementarity_scaling);
   }
} // namespace
