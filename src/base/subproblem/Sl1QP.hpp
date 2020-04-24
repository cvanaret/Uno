#ifndef Sl1QP_H
#define Sl1QP_H

#include "Subproblem.hpp"
#include "QPSolver.hpp"

/*! \class QPApproximation
 * \brief QP local approximation
 *
 *  Quadratic approximation
 */
class Sl1QP : public Subproblem {
public:
    /*!
     *  Constructor
     * 
     * \param solver: solver that solves the subproblem
     */
    Sl1QP(QPSolver& solver);

    Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, bool use_trust_region) override;

    SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double trust_region_radius = INFINITY) override;
    SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius = INFINITY) override;
    
    void compute_optimality_measures(Problem& problem, Iterate& iterate) override;
    void compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution);
    
    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) override;
    bool phase_1_required(SubproblemSolution& solution) override;

    /* use a reference to allow polymorphism */
    QPSolver& solver; /*!< Subproblem solver */

private:
    /* problem reformulation */
    // equality constraints c(x) = 0 are reformulated as c(x) = p - n
    // equality constraints l <= c(x) = u are reformulated as c(x) - s = p - n
    std::map<int, int> slack_variables; // s
    std::map<int, int> positive_part_variables; // p
    std::map<int, int> negative_part_variables; // n

    std::vector<Range> generate_variables_bounds(Problem& problem, Iterate& current_iterate, double trust_region_radius);
    void recover_active_set(SubproblemSolution& solution, std::vector<Range>& variables_bounds);
};

#endif // Sl1QP_H
