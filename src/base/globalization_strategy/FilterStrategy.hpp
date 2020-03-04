#ifndef FILTERSTRATEGY_H
#define FILTERSTRATEGY_H

#include <iostream>
#include <memory>
#include "TwoPhaseStrategy.hpp"
#include "Filter.hpp"

/*! \class FilterStrategy
 * \brief Step acceptance strategy based on a filter
 *
 *  Strategy that accepts or declines a trial step
 */
class FilterStrategy : public TwoPhaseStrategy {
public:
    /*!
     *  Constructor that takes an optimization problem, filters for restoration and optimality, and a set of constants
     */
    FilterStrategy(Subproblem& subproblem, std::shared_ptr<Filter> filter_optimality,
        std::shared_ptr<Filter> filter_restoration, TwoPhaseConstants& constants, double tolerance);

    /* use pointers to allow polymorphism */
    std::shared_ptr<Filter> filter_optimality; /*!< Filter for the optimality phase */
    std::shared_ptr<Filter> filter_restoration; /*!< Filter for the restoration phase */

    /*!
     *  Check the validity of a step
     *  Implements the purely virtual method of the superclass
     */
    bool check_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length = 1.) override;
    Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool use_trust_region) override;
    void compute_measures(Problem& problem, Iterate& iterate) override;

private:
    void switch_phase(Problem& problem, SubproblemSolution& solution, Iterate& current_iterate, Iterate& trial_iterate);
    OptimalityStatus compute_status(Problem& problem, Iterate& current_iterate, double step_norm);
};

#endif // FILTERSTRATEGY_H
