// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRUSTREGIONSTRATEGY_H
#define UNO_TRUSTREGIONSTRATEGY_H

#include "GlobalizationMechanism.hpp"

class TrustRegionStrategy : public GlobalizationMechanism {
public:
   TrustRegionStrategy(Statistics& statistics, ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);

   void initialize(Iterate& initial_iterate) override;
   Iterate compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate) override;

private:
   double radius; /*!< Current trust region radius */
   const double increase_factor;
   const double decrease_factor;
   const double aggressive_decrease_factor;
   const double activity_tolerance;
   const double minimum_radius;
   const double radius_reset_threshold;

   Iterate assemble_trial_iterate(const Model& model, Iterate& current_iterate, const Direction& direction);
   void possibly_increase_radius(double step_norm);
   void decrease_radius(double step_norm);
   void decrease_radius();
   void decrease_radius_aggressively();
   void reset_radius();
   void reset_active_trust_region_multipliers(const Model& model, const Direction& direction, Iterate& trial_iterate) const;
   void set_statistics(Statistics& statistics, const Direction& direction);
   [[nodiscard]] bool termination() const;
   void print_iteration();
};

#endif // UNO_TRUSTREGIONSTRATEGY_H
