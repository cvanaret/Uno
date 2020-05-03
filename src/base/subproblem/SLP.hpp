#ifndef SLP_H
#define SLP_H

#include "ActiveSetMethod.hpp"

/*! \class SLP
 * \brief LP local approximation
 *
 *  Linear approximation
 */
class SLP : public ActiveSetMethod {
public:
    /*!
     *  Constructor
     * 
     * \param solver: solver that solves the subproblem
     */
    SLP(QPSolver& solver);

    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution);
    bool phase_1_required(SubproblemSolution& solution);

private:
    void evaluate_optimality_iterate(Problem& problem, Iterate& current_iterate);
    // phase 1
    void evaluate_feasibility_iterate(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution);
    // call subproblem solver
    SubproblemSolution solve_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, Iterate& current_iterate, std::vector<double>& d0);
};

#endif // SLP_H
