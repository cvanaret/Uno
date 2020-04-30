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
        PenaltyMeritFunction(Subproblem& subproblem, double tolerance);

        /*!
         *  Check the validity of a step
         */
        bool check_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) override;
        Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool use_trust_region) override;

    private:
        double eta;
        
        double compute_linear_model(Problem& problem, SubproblemSolution& solution);
        //double compute_error(Problem& problem, Iterate& current_iterate, Multipliers& multipliers, double penalty_parameter);
        OptimalityStatus compute_status(Problem& problem, Iterate& current_iterate, double step_norm, double objective_multiplier);
};

#endif // PENALTYMERITFUNCTION_H
