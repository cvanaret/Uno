// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BACKTRACKINGLINESEARCH_H
#define UNO_BACKTRACKINGLINESEARCH_H

#include "GlobalizationMechanism.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class Direction;
   class PredictedReductionModels;
   class ProgressMeasures;

   class BacktrackingLineSearch : public GlobalizationMechanism {
   public:
      BacktrackingLineSearch(const Model& model, Options& options);
      ~BacktrackingLineSearch() override = default;

      void initialize(Statistics& statistics, const Model& model, Iterate& current_iterate, EvaluationCache& evaluation_cache,
         Options& options) override;
      void compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         EvaluationCache& evaluation_cache, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;

      [[nodiscard]] std::string get_name() const override;

   private:
      const double backtracking_ratio;
      const bool scale_duals_with_step_length;
      // minimum step length (Eq. 23 in IPOPT paper)
      const double delta;
      const double gamma_alpha;
      const double gamma_theta;
      const double gamma_phi;
      const double s_theta;
      const double s_phi;
      const double theta_min;
      // tiny directions
      size_t number_consecutive_tiny_directions{0};
      const size_t consecutive_tiny_directions_threshold{2}; // TODO add option
      // second-order corrections
      const size_t SOC_max_iterations;
      const double SOC_infeasibility_fraction;

      void assemble_trial_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) const;
      [[nodiscard]] double compute_minimum_step_length(const PredictedReductionModels& predicted_reduction_models,
         const Iterate& current_iterate, double objective_multiplier) const;
      [[nodiscard]] bool backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, EvaluationCache& evaluation_cache,
         const PredictedReductionModels& predicted_reduction_models, double minimum_step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks);
      [[nodiscard]] bool is_tiny_direction(const Iterate& current_iterate, const Direction& direction) const;
      [[nodiscard]] bool compute_second_order_directions(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, EvaluationCache& evaluation_cache,
         const ProgressMeasures& predicted_reductions, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) const;
      static void check_unboundedness(const Direction& direction);
   };
} // namespace

#endif // UNO_BACKTRACKINGLINESEARCH_H