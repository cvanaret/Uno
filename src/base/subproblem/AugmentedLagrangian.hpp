#ifndef AUGMENTEDLAGRANGIAN_H
#define AUGMENTEDLAGRANGIAN_H

#include <map>
#include "Subproblem.hpp"
#include "LBFGSB.hpp"

/*! \class AugmentedLagrangian
 * \brief Augmented Lagrangian subproblem method
 *
 *  Strategy that computes a (descend) direction by approx. minimization of the augmented Lagrangian
 */
class AugmentedLagrangian : public Subproblem {
    public:
        /*!
         *  Constructor
         * 
         * \param solver: solver that computes the step
         */
        AugmentedLagrangian();

        Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, int number_constraints, bool use_trust_region) override;

        SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius) override;
        SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, SubproblemSolution& phase_II_solution) override;
        SubproblemSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) override;
        void compute_measures(Problem& problem, Iterate& iterate) override;
        bool phase_1_required(SubproblemSolution& solution) override;

        LBFGSB solver; /*!< Solver that computes the step */
};

#endif // AUGMENTEDLAGRANGIAN_H
