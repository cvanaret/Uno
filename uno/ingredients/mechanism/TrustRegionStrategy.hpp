#ifndef TRUSTREGION_H
#define TRUSTREGION_H

#include "GlobalizationMechanism.hpp"

/*! \class TrustRegion
 * \brief Trust region
 *
 *  Trust region strategy
 */
class TrustRegionStrategy : public GlobalizationMechanism {
public:
   TrustRegionStrategy(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);

   void initialize(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& first_iterate) override;
   std::tuple<Iterate, double> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, const Scaling& scaling,
         Iterate& current_iterate) override;

private:
   double radius; /*!< Current trust region radius */
   const double increase_factor;
   const double decrease_factor;
   const double activity_tolerance;
   const double min_radius;

   static void rectify_active_set(Direction& direction, double radius);
   void add_statistics(Statistics& statistics, const Direction& direction);
   static void check_unboundedness(const Direction& direction);
   bool termination();
   void print_iteration();
};

#endif // TRUSTREGION_H
