#ifndef GLOBALIZATIONMECHANISM_H
#define GLOBALIZATIONMECHANISM_H

#include "Problem.hpp"
#include "GlobalizationStrategy.hpp"

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
        GlobalizationMechanism(GlobalizationStrategy& globalization_strategy, int max_iterations);
        virtual ~GlobalizationMechanism();

        virtual Iterate compute_iterate(Problem& problem, Iterate& current_iterate) = 0;
        virtual Iterate initialize(Problem& problem, std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers) = 0;

        /* references to allow polymorphism */
        GlobalizationStrategy& globalization_strategy; /*!< Strategy to accept or reject a step */
        int max_iterations; /*!< Maximum number of iterations */
        int number_iterations; /*!< Current number of iterations */
};

#endif // GLOBALIZATIONMECHANISM_H
