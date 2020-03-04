#ifndef SQP_H
#define SQP_H

#include "Subproblem.hpp"
#include "QPSolver.hpp"

/*! \class QPApproximation
 * \brief QP local approximation
 *
 *  Quadratic approximation
 */
class SQP : public Subproblem {
public:
    /*!
     *  Constructor
     * 
     * \param solver: solver that solves the subproblem
     */
    SQP(QPSolver& solver);

    Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, int number_constraints, bool use_trust_region) override;

    SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) override;
    SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) override;
    SubproblemSolution compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, double penalty_parameter, PenaltyDimensions penalty_dimensions) override;
    void compute_measures(Problem& problem, Iterate& iterate) override;
    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) override;
    bool phase_1_required(SubproblemSolution& solution) override;

    /* use a reference to allow polymorphism */
    QPSolver& solver; /*!< Solver that solves the subproblem */

private:
    // phase 1
    std::vector<double> generate_feasibility_multipliers(Problem& problem, std::vector<double>& current_constraint_multipliers, ConstraintPartition& constraint_partition);
    std::vector<Range> generate_feasibility_bounds(Problem& problem, std::vector<double>& current_constraints, ConstraintPartition& constraint_partition);
    void set_feasibility_objective_(Problem& problem, Iterate& current_iterate, ConstraintPartition& constraint_partition);
};

#endif // SQP_H
