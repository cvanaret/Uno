#ifndef TRUSTLINESEARCH_H
#define TRUSTLINESEARCH_H

#include "GlobalizationMechanism.hpp"

/*! \class LineSearch
 * \brief Line-search
 *
 *  Line-search strategy
 */
class TrustLineSearch : public GlobalizationMechanism {
    public:
        /*!
         *  Constructor
         */
        TrustLineSearch(GlobalizationStrategy& globalization_strategy, double initial_radius, int max_iterations = 30, double ratio = 0.5);

        /*!
         *  Compute the next iterate from a given point
         * 
         * \param problem: optimization problem
         * \param current_iterate: current point and its evaluations
         */
        Iterate compute_acceptable_iterate(Problem& problem, Iterate& current_iterate);
        Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;

        /* ratio of step length update in ]0, 1[ */
        double ratio;
        double radius; /*!< Current trust region radius */

    private:
        bool termination(bool is_accepted, int iteration);
        std::vector<Range> compute_variables_bounds(Problem& problem, Iterate& current_iterate, double radius);
        void correct_multipliers(Problem& problem, SubproblemSolution& solution);

        double activity_tolerance_;
};

#endif // TRUSTLINESEARCH_H
