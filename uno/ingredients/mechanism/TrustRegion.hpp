#ifndef TRUSTREGION_H
#define TRUSTREGION_H

#include "GlobalizationMechanism.hpp"

/*! \class TrustRegion
 * \brief Trust region
 *
 *  Trust region strategy
 */
class TrustRegion : public GlobalizationMechanism {
public:
   TrustRegion(ConstraintRelaxationStrategy& constraint_relaxation_strategy, double initial_radius, int max_iterations);

   void initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) override;
   std::tuple<Iterate, double> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;

private:
   double radius; /*!< Current trust region radius */
   const double increase_factor{2.};
   const double decrease_factor{2.};
   const double activity_tolerance{1e-6};
   const double min_radius{1e-16};
   const double full_step_length{1.};

   static void rectify_active_set(Direction& direction, double radius);
   void add_statistics(Statistics& statistics, const Direction& direction);
   static void check_unboundedness(const Direction& direction);
   bool termination();
   void print_iteration();
};

#endif // TRUSTREGION_H
