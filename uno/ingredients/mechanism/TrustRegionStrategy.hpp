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

   void initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) override;
   std::tuple<Iterate, double> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;

private:
   double radius; /*!< Current trust region radius */
   const double increase_factor;
   const double decrease_factor;
   const double activity_tolerance;
   const double min_radius;

   static void rectify_active_set(Direction& direction, double radius);
   void add_statistics(Statistics& statistics, const Direction& direction);
   bool termination() const;
   void print_iteration();
};

#endif // UNO_TRUSTREGIONSTRATEGY_H
