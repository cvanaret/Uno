#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <vector>
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Phase.hpp"
#include "SubproblemSolution.hpp"
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
        Subproblem();
        virtual ~Subproblem();

        virtual SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) = 0;
        virtual SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) = 0;
        virtual SubproblemSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, double penalty_parameter, PenaltyDimensions penalty_dimensions) = 0;

        virtual Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, int number_constraints, bool use_trust_region) = 0;
        virtual void compute_measures(Problem& problem, Iterate& iterate) = 0;
        virtual double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) = 0;
        virtual bool phase_1_required(SubproblemSolution& solution) = 0;
        
        static std::vector<Range> generate_constraints_bounds(Problem& problem, std::vector<double>& current_constraints);
        
        int number_subproblems_solved;
};

#endif // SUBPROBLEM_H
