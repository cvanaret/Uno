#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <vector>
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Phase.hpp"
#include "LocalSolution.hpp"
#include "Constraint.hpp"

/*! \class Subproblem
 * \brief Subproblem
 *
 *  Local approximation of a nonlinear optimization problem (virtual class) 
 */
class Subproblem {
    public:
        /*!
         *  Constructor
         * 
         * \param solver: solver that solves the subproblem
         * \param name: name of the strategy
         */
        Subproblem(std::string name);
        virtual ~Subproblem();

        virtual LocalSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius) = 0;

        virtual LocalSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, LocalSolution& phase_II_solution) = 0;

        virtual LocalSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) = 0;

        virtual void initialize(Problem& problem, Iterate& current_iterate, int number_variables, int number_constraints, bool use_trust_region) = 0;

        std::string name; /*!< Name of the strategy */
        int number_subproblems_solved;
};

#endif // SUBPROBLEM_H
