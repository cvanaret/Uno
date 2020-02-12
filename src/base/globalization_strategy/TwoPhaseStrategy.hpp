#ifndef TWOPHASESTRATEGY_H
#define TWOPHASESTRATEGY_H

#include "GlobalizationStrategy.hpp"

/*! \class TwoPhaseConstants
 * \brief Constants for filter and tube strategies
 *
 *  Set of constants to control the filter and tube strategies
 */
struct TwoPhaseConstants {
    double Sigma; /*!< Sufficient reduction constant */
    double Delta; /*!< Switching constant */
    double ubd;
    double fact;
};

/*! \class GlobalizationStrategy
 * \brief Step acceptance strategy
 *
 *  Strategy that accepts or declines a trial step (virtual class)
 */
class TwoPhaseStrategy : public GlobalizationStrategy {
    public:
        TwoPhaseStrategy(Subproblem& subproblem, TwoPhaseConstants& constants, double tolerance);

        SubproblemSolution compute_step(Problem& problem, Iterate& current_iterate, double radius) override;

        Phase phase; /*!< Current phase (optimality or feasibility restoration) */
        TwoPhaseConstants constants; /*!< Set of constants */

    protected:
        void update_restoration_multipliers(Iterate& trial_iterate, ConstraintPartition& constraint_partition);
};

#endif // TWOPHASESTRATEGY_H
