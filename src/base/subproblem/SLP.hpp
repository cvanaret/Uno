#ifndef SLP_H
#define SLP_H

#include "Subproblem.hpp"
#include "QPSolver.hpp"

/*! \class SLP
 * \brief LP local approximation
 *
 *  Linear approximation
 */
class SLP : public Subproblem {
public:
    /*!
     *  Constructor
     * 
     * \param solver: solver that solves the subproblem
     */
    SLP(LPSolver& solver);

    Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, int number_constraints, std::vector<Range>& variables_bounds, bool use_trust_region) override;

    SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) override;
    SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) override;
    //SubproblemSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, double penalty_parameter, PenaltyDimensions penalty_dimensions) override;
    void compute_measures(Problem& problem, Iterate& iterate) override;
    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) override;
    bool phase_1_required(SubproblemSolution& solution) override;

    /* use a reference to allow polymorphism */
    LPSolver& solver; /*!< Solver that solves the subproblem */

private:
    void set_infeasibility_objective_(Problem& problem, Iterate& current_iterate, ConstraintPartition& constraint_partition);
};

#endif // SLP_H
