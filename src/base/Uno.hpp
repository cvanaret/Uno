#ifndef UNO_H
#define UNO_H

#include "Problem.hpp"
#include "GlobalizationMechanism.hpp"

struct Result {
    Iterate solution;
    int number_variables;
    int number_constraints;
    int iteration;
    double cpu_time;
    int objective_evaluations;
    int constraint_evaluations;
    int jacobian_evaluations;
    int hessian_evaluations;
    int number_subproblems_solved;

    void display(bool print_solution);
};

/*! \class Argonot
 * \brief Argonot
 *
 *  UNO solver
 */
class Uno {
public:
    /*!
     *  Constructor
     * 
     * \param globalization_strategy: strategy to promote global convergence
     * \param tolerance: tolerance for termination criteria
     */
    Uno(GlobalizationMechanism& globalization_mechanism, int max_iterations);

    /*!
     *  Solve a given problem with initial primal and dual variables
     * 
     * \param problem: optimization problem
     * \param x: primal variables
     * \param multipliers: Lagrange multipliers/dual variables
     */
    Result solve(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool preprocessing);
    void preprocessing(Problem& problem, std::vector<double>& x, Multipliers& multipliers);

    GlobalizationMechanism& globalization_mechanism; /*!< Step control strategy (trust region or line-search) */
    int max_iterations; /*!< Maximum number of iterations */

private:
    bool termination_criterion_(TerminationStatus is_optimal, int iteration);
    TerminationStatus optimality_test_(Problem& problem, Phase& phase, Iterate& current_iterate);
};

#endif // UNO_H
