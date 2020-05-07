#ifndef FILTERSTRATEGY_H
#define FILTERSTRATEGY_H

#include <iostream>
#include <memory>
#include "Filter.hpp"

/*! \class TwoPhaseConstants
 * \brief Constants for filter and tube strategies
 *
 *  Set of constants to control the filter and tube strategies
 */
struct FilterStrategyParameters {
    double Sigma; /*!< Sufficient reduction constant */
    double Delta; /*!< Switching constant */
    double ubd;
    double fact;
};

/*! \class FilterStrategy
 * \brief Step acceptance strategy based on a filter
 *
 *  Strategy that accepts or declines a trial step
 */
class FilterStrategy : public GlobalizationStrategy {
public:
    /*!
     *  Constructor that takes an optimization problem, filters for restoration and optimality, and a set of constants
     */
    FilterStrategy(Subproblem& subproblem, std::shared_ptr<Filter> filter_optimality, std::shared_ptr<Filter> filter_restoration, FilterStrategyParameters& strategy_parameters, double tolerance);

    /* use pointers to allow polymorphism */
    std::shared_ptr<Filter> filter_optimality; /*!< Filter for the optimality phase */
    std::shared_ptr<Filter> filter_restoration; /*!< Filter for the restoration phase */

    /*!
     *  Check the validity of a step
     *  Implements the purely virtual method of the superclass
     */
    bool check_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length = 1.) override;
    Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool use_trust_region) override;
    
private:
    Phase current_phase; /*!< Current phase (optimality or feasibility restoration) */
    FilterStrategyParameters parameters; /*!< Set of constants */
    
    void switch_phase(Problem& problem, SubproblemSolution& solution, Iterate& current_iterate, Iterate& trial_iterate);
    void update_restoration_multipliers(Iterate& trial_iterate, ConstraintPartition& constraint_partition);
};

#endif // FILTERSTRATEGY_H
