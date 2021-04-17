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
    /*!
     *  Constructor
     * 
     * \param direction_computation: strategy to compute a descent direction
     * \param step_accept: strategy to accept or reject a step
     * \param initial_radius: initial trust region radius
     */
    TrustRegion(GlobalizationStrategy& globalization_strategy, double tolerance, double initial_radius, int max_iterations);

    Iterate initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;
    Iterate compute_acceptable_iterate(Statistics& statistics, Problem& problem, Iterate& current_iterate) override;

    static void correct_active_set(Direction& direction, double radius);

    double radius; /*!< Current trust region radius */

private:
    double activity_tolerance_;

    bool termination_(bool is_accepted);
    void print_iteration_();
    void print_acceptance_() override;
    void print_warning_(const char* message);
};

#endif // TRUSTREGION_H
