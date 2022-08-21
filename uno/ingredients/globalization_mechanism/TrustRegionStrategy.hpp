// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRUSTREGIONSTRATEGY_H
#define UNO_TRUSTREGIONSTRATEGY_H

#include "GlobalizationMechanism.hpp"

/*! \class TrustRegionStrategy
 * \brief Trust region strategy
 *
 *  Trust region strategy
 */
class TrustRegionStrategy : public GlobalizationMechanism {
public:
   TrustRegionStrategy(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);

   void initialize(Statistics& statistics, Iterate& first_iterate) override;
   std::tuple<Iterate, double> compute_acceptable_iterate(Statistics& statistics, Iterate& current_iterate) override;

private:
   double radius; /*!< Current trust region radius */
   const double increase_factor;
   const double decrease_factor;
   const double activity_tolerance;
   const double min_radius;
   // statistics table
   int statistics_TR_radius_column_order;

   void increase_radius(double step_norm);
   void decrease_radius(double step_norm);
   static void rectify_multipliers(Direction& direction, double radius);
   void set_statistics(Statistics& statistics, const Direction& direction);
   [[nodiscard]] bool termination() const;
   void print_iteration();
};

#endif // UNO_TRUSTREGIONSTRATEGY_H
