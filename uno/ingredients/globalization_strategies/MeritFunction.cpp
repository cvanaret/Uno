// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MeritFunction.hpp"
#include "ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/Symbols.hpp"

namespace uno {
   MeritFunction::MeritFunction(const Options& options):
         GlobalizationStrategy(options),
         sufficient_infeasibility_decrease_ratio(options.get_double("sufficient_infeasibility_decrease_ratio")),
         infeasibility_reduction_fraction(options.get_double("merit_infeasibility_reduction_fraction")) {
   }

   void MeritFunction::initialize(Statistics& statistics, const Iterate& initial_iterate) {
      this->smallest_known_infeasibility = initial_iterate.progress.infeasibility;
      this->inverse_penalty_parameter = 1.;
      statistics.add_column("Penalty", Statistics::double_width, 2, /* is_extended = */ true);
   }

   Decision MeritFunction::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reductions, double objective_multiplier) {
      // the inverse penalty parameter is local to the merit function: the problem (objective multiplier) is unchanged
      this->update_inverse_penalty_parameter(predicted_reductions, objective_multiplier);
      statistics.set("Penalty", this->inverse_penalty_parameter);

      // predicted reduction with all contributions. This quantity should be positive (= negative directional derivative)
      const double constrained_predicted_reduction = this->merit_function(predicted_reductions, objective_multiplier);
      DEBUG << "Constrained predicted reduction: " << constrained_predicted_reduction << '\n';
      if (constrained_predicted_reduction <= 0.) {
         DEBUG  << "The direction is not a descent direction for the merit function\n";
      }
      // compute current exact penalty
      const double current_merit_value = this->merit_function(current_progress, objective_multiplier);
      const double trial_merit_value = this->merit_function(trial_progress, objective_multiplier);
      const double actual_reduction = this->compute_actual_reduction(current_merit_value, trial_merit_value);
      DEBUG << "Current merit: " << current_merit_value << '\n';
      DEBUG << "Trial merit:   " << trial_merit_value << '\n';
      DEBUG << "Actual reduction: " << current_merit_value << " - " << trial_merit_value << " = " << actual_reduction << '\n';

      // Armijo sufficient decrease condition
      if (this->armijo_sufficient_decrease(constrained_predicted_reduction, actual_reduction)) {
         DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
         this->smallest_known_infeasibility = std::min(this->smallest_known_infeasibility, trial_progress.infeasibility);
         statistics.set("Status", std::string(symbols::check) + " (Armijo)");
         return Decision::accept();
      }
      DEBUG << "Trial iterate was rejected by Armijo condition\n";
      statistics.set("Status", std::string(symbols::fail) + " (Armijo)");
      return Decision::reject(RejectionReason::SUFFICIENT_DECREASE);
   }

   bool MeritFunction::is_infeasibility_sufficiently_reduced(const Iterate& trial_iterate, double /*reference_infeasibility*/) const {
      // if the trial infeasibility sufficiently decreases the best known infeasibility
      return (trial_iterate.progress.infeasibility <= this->sufficient_infeasibility_decrease_ratio * this->smallest_known_infeasibility);
   }

   void MeritFunction::reset() {
      this->inverse_penalty_parameter = 1.;
   }

   void MeritFunction::avoid_cycling_back_to(const ProgressMeasures& /*current_progress*/) {
   }

   std::string MeritFunction::get_name() const {
      return "merit";
   }

   // protected member functions

   // merit function τ (objective + auxiliary) + infeasibility, i.e. the exact penalty (objective + auxiliary) + 1/τ infeasibility
   // scaled by τ. The whole optimality measure is scaled: the direction is computed for objective + auxiliary
   double MeritFunction::merit_function(const ProgressMeasures& progress, double objective_multiplier) const {
      return this->inverse_penalty_parameter * (progress.objective(objective_multiplier) + progress.auxiliary) +
         progress.infeasibility;
   }

   // decrease τ such that the predicted reduction of the merit function satisfies pred(τ) ≥ σ Δinfeasibility, σ ∈ (0, 1)
   // (Nocedal & Wright, eq. (18.36) with σ = 0, in inverse form). The condition is linear in τ:
   // pred(τ) = τ Δoptimality + Δinfeasibility ≥ σ Δinfeasibility <=> τ (-Δoptimality) ≤ (1 - σ) Δinfeasibility
   void MeritFunction::update_inverse_penalty_parameter(const ProgressMeasures& predicted_reductions, double objective_multiplier) {
      const double predicted_optimality_reduction = predicted_reductions.objective(objective_multiplier) + predicted_reductions.auxiliary;
      const double predicted_infeasibility_reduction = predicted_reductions.infeasibility;
      // if the infeasibility does not decrease (Δinfeasibility ≤ 0), decreasing τ cannot help
      // if the optimality measure decreases (Δoptimality ≥ 0), the condition holds for all τ ≥ 0
      if (0. < predicted_infeasibility_reduction && predicted_optimality_reduction < 0.) {
         const double threshold = (1. - this->infeasibility_reduction_fraction) * predicted_infeasibility_reduction /
            -predicted_optimality_reduction;
         if (threshold < this->inverse_penalty_parameter) {
            this->inverse_penalty_parameter = threshold;
            DEBUG << "Inverse penalty parameter decreased to " << this->inverse_penalty_parameter << '\n';
         }
      }
   }
} // namespace