// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BACKTRACKINGLINESEARCH_H
#define UNO_BACKTRACKINGLINESEARCH_H

#include <optional>
#include "GlobalizationMechanism.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   // forward declaration
   class PredictedReductionModels;
   class ProgressMeasures;

   class BacktrackingLineSearch : public GlobalizationMechanism {
   public:
      BacktrackingLineSearch(const Model& model, const Options& options, std::vector<OptionOverride>& option_overrides);
      ~BacktrackingLineSearch() override = default;

      void initialize(Statistics& statistics, const Model& model, Iterate& current_iterate, EvaluationCache& evaluation_cache,
         const Options& options, std::vector<OptionOverride>& option_overrides) override;
      void compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         EvaluationCache& evaluation_cache, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;
      void postprocess_solution(Iterate& iterate, Evaluations& evaluations) const override;

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
      // IPOPT hard-codes this absolute threshold in BacktrackingLineSearch::DetectTinyStep (independent of theta_min)
      static constexpr double tiny_step_infeasibility_threshold{1e-4};
      // second-order corrections
      const size_t SOC_max_iterations;
      const double SOC_infeasibility_fraction;
      // tiny directions (IPOPT): tolerance (multiple of the machine epsilon) on the primal direction, tolerance on the
      // constraint dual direction, and whether the previous iteration was a tiny step with a small dual direction
      const double tiny_step_factor;
      static constexpr double tiny_step_dual_tolerance{1e-2}; // IPOPT's tiny_step_y_tol
      const double primal_tolerance;
      bool previous_tiny_step_with_small_dual_direction{false};
      const bool print_minor_iterations;

      // watchdog (Chamberlain et al., 1982; IPOPT's BacktrackingLineSearch): after watchdog_shortened_iterations_trigger
      // consecutive iterations with a shortened step, the current iterate becomes the watchdog (reference) iterate. Then,
      // for at most watchdog_max_trial_iterations iterations, the full step is accepted without globalization test. As
      // soon as a trial iterate is acceptable wrt the watchdog iterate, the watchdog stops. Otherwise, the algorithm
      // reverts to the watchdog iterate and backtracks along the watchdog direction. 0 disables the watchdog
      const size_t watchdog_shortened_iterations_trigger;
      const size_t watchdog_max_trial_iterations;
      size_t number_shortened_iterations{0};
      bool in_watchdog{false};
      size_t watchdog_trial_iterations{0};
      std::optional<Iterate> watchdog_iterate{};
      std::optional<Evaluations> watchdog_evaluations{};
      Direction watchdog_direction{};
      size_t last_number_parameterization_updates{0};
      bool last_step_was_shortened{false};

      void assemble_trial_iterate(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) const;
      [[nodiscard]] double compute_minimum_step_length(const PredictedReductionModels& predicted_reduction_models,
         const Iterate& current_iterate) const;
      [[nodiscard]] bool backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, EvaluationCache& evaluation_cache,
         const PredictedReductionModels& predicted_reduction_models, double minimum_step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks, bool skip_first_trial = false);
      void update_watchdog_state(const Iterate& current_iterate, const Direction& direction, const Evaluations& current_evaluations);
      [[nodiscard]] bool take_watchdog_step(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, EvaluationCache& evaluation_cache,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks);
      void stop_watchdog();
      [[nodiscard]] bool try_tiny_step(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, EvaluationCache& evaluation_cache, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks);
      [[nodiscard]] bool is_tiny_direction(const Iterate& current_iterate, const Direction& direction) const;
      [[nodiscard]] bool compute_second_order_directions(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, EvaluationCache& evaluation_cache,
         const ProgressMeasures& predicted_reductions, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) const;
      static void check_unboundedness(const Direction& direction);
   };
} // namespace

#endif // UNO_BACKTRACKINGLINESEARCH_H