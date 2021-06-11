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
   double radius; /*!< Current trust region radius */

   /*!
    *  Constructor
    *
    * \param direction_computation: strategy to compute a descent direction
    * \param step_accept: strategy to accept or reject a step
    * \param initial_radius: initial trust region radius
    */
   TrustRegion(ConstraintRelaxationStrategy& constraint_relaxation_strategy, double initial_radius, int max_iterations);

   Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;
   std::pair<Iterate, Direction> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   static void rectify_active_set(Direction& direction, double radius);

private:
   double activity_tolerance_;

   void add_statistics(Statistics& statistics, const Direction& direction);
   bool termination();
   void print_iteration();
};

#endif // TRUSTREGION_H
