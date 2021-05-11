#ifndef GLOBALIZATIONMECHANISM_H
#define GLOBALIZATIONMECHANISM_H

#include "Problem.hpp"
#include "GlobalizationStrategy.hpp"
#include "Statistics.hpp"

/*! \class GlobalizationMechanism
 * \brief Step control strategy
 *
 *  Strategy that promotes global convergence
 */
class GlobalizationMechanism {
public:
    /*!
     *  Constructor
     * 
     * \param direction_computation: strategy to compute a descent direction
     * \param step_accept: strategy to accept or reject a step
     */
    GlobalizationMechanism(GlobalizationStrategy& globalization_strategy, double tolerance, int max_iterations);
    virtual ~GlobalizationMechanism();

    virtual Iterate initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) = 0;
    virtual Iterate compute_acceptable_iterate(Statistics& statistics, Problem& problem, Iterate& current_iterate) = 0;
    std::optional<std::pair<Iterate, Direction> > find_first_acceptable_direction_(Statistics& statistics, Problem& problem, Iterate& current_iterate, std::vector<Direction>& directions, double step_length);

    /* references to allow polymorphism */
    GlobalizationStrategy& globalization_strategy; /*!< Strategy to accept or reject a step */
    double tolerance; /*!< Tolerance of the termination criteria */
    int max_iterations; /*!< Maximum number of iterations */
    int number_iterations; /*!< Current number of iterations */

protected:
    TerminationStatus compute_termination_status_(Problem& problem, Iterate& current_iterate, double step_norm, double objective_multiplier);
    virtual void print_acceptance_() = 0;
};

#endif // GLOBALIZATIONMECHANISM_H
