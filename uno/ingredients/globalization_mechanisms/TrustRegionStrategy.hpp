// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRUSTREGIONSTRATEGY_H
#define UNO_TRUSTREGIONSTRATEGY_H

#include "GlobalizationMechanism.hpp"

namespace uno {
   class TrustRegionStrategy : public GlobalizationMechanism {
   public:
      TrustRegionStrategy(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);

      void initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) override;
      void compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
            WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;

   private:
      double radius; /*!< Current trust region radius */
      const double increase_factor;
      const double decrease_factor;
      const double aggressive_decrease_factor;
      const double activity_tolerance;
      const double minimum_radius;
      const double radius_reset_threshold;
      const double tolerance;

      bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
            WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks);
      void possibly_increase_radius(double step_norm);
      void decrease_radius(double step_norm);
      void decrease_radius();
      void decrease_radius_aggressively();
      void reset_radius();
      void reset_active_trust_region_multipliers(const Model& model, const Direction& direction, Iterate& trial_iterate) const;
      bool check_termination_with_small_step(Iterate& trial_iterate) const;
      void set_trust_region_statistics(Statistics& statistics, size_t number_iterations) const;
      void set_statistics(Statistics& statistics, const Direction& direction) const;
      void set_statistics(Statistics& statistics, const Iterate& trial_iterate, const Direction& direction) const;
   };
} // namespace

#endif // UNO_TRUSTREGIONSTRATEGY_H
