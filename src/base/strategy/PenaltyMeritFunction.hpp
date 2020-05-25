#ifndef PENALTYMERITFUNCTION_H
#define PENALTYMERITFUNCTION_H

#include <vector>
#include "GlobalizationStrategy.hpp"
#include "Constraint.hpp"

/*! \class PenaltyStrategy
 * \brief Step acceptance strategy based on a penalty method
 *
 *  Strategy that accepts or declines a trial step
 */
class PenaltyMeritFunction : public GlobalizationStrategy {
    public:
        /*!
         *  Constructor that takes an optimization problem and a set of constants
         */
        PenaltyMeritFunction(Subproblem& subproblem);

        /*!
         *  Check the validity of a step
         */
        bool check_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) override;
        Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;

    private:
        double decrease_fraction_;
        
        double compute_linear_model_(Problem& problem, SubproblemSolution& solution);
};

#endif // PENALTYMERITFUNCTION_H
