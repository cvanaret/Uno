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
    TrustRegion(GlobalizationStrategy& globalization_strategy, double initial_radius, int max_iterations = 100);

    Iterate compute_iterate(Problem& problem, Iterate& current_iterate) override;
    Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;

    double radius; /*!< Current trust region radius */

private:
    double activity_tolerance_;

    std::vector<Range> compute_trust_region(Iterate& current_iterate, double radius);
    void correct_multipliers(Problem& problem, SubproblemSolution& solution);
    bool termination(bool is_accepted);
    void print_iteration();
    void print_acceptance(double solution_norm);
    void print_warning(const char* message);
};

#endif // TRUSTREGION_H
