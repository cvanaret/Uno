#ifndef FILTERSTRATEGY_H
#define FILTERSTRATEGY_H

#include <iostream>
#include <memory>
#include "GlobalizationStrategy.hpp"
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
    FilterStrategy(Subproblem& subproblem, FilterStrategyParameters strategy_constants, std::map<std::string, std::string> options);

    /* use pointers to allow polymorphism */
    std::shared_ptr<Filter> filter_optimality; /*!< Filter for the optimality phase */
    std::shared_ptr<Filter> filter_restoration; /*!< Filter for the restoration phase */

    Iterate initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;
    std::optional<Iterate> check_acceptance(Statistics& statistics, Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) override;
    
private:
    Phase current_phase_; /*!< Current phase (optimality or feasibility restoration) */
    FilterStrategyParameters parameters_; /*!< Set of constants */
    
    void switch_phase_(Problem& problem, Direction& direction, Iterate& current_iterate, Iterate& trial_iterate);
    void update_restoration_multipliers_(Iterate& trial_iterate, ConstraintPartition& constraint_partition);
};

#endif // FILTERSTRATEGY_H
