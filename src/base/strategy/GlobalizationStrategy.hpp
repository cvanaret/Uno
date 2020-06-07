#ifndef GLOBALIZATIONSTRATEGY_H
#define GLOBALIZATIONSTRATEGY_H

#include <cmath>
#include "Problem.hpp"
#include "Subproblem.hpp"
#include "Iterate.hpp"
#include "Direction.hpp"
#include "Constraint.hpp"

/*! \class GlobalizationStrategy
 * \brief Step acceptance strategy
 *
 *  Strategy that accepts or declines a trial step (virtual class)
 */
class GlobalizationStrategy {
public:
    /*!
     *  Constructor that takes an optimization problem and a tolerance
     * 
     * \param problem: optimization problem
     * \param constants: set of constants
     */
    GlobalizationStrategy(Subproblem& subproblem);
    virtual ~GlobalizationStrategy();

    Subproblem& subproblem;

    /*!
     *  Check the validity of a step
     *  Purely virtual method (only implemented in subclasses)
     */
    virtual bool check_step(Problem& problem, Iterate& current_iterate, Direction& direction, double step_length = 1.) = 0;
    virtual Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers) = 0;
};

#endif // GLOBALIZATIONSTRATEGY_H
