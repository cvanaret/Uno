#ifndef ACTIVESETMETHOD_H
#define ACTIVESETMETHOD_H

#include <memory>
#include "Subproblem.hpp"
#include "QPSolver.hpp"
#include "HessianEvaluation.hpp"

class ActiveSetMethod : public Subproblem {
public:
    ActiveSetMethod(Problem& problem, std::shared_ptr<QPSolver> solver, bool scale_residuals);

    virtual Iterate evaluate_initial_point(Problem& problem, std::vector<double>& x, Multipliers& multipliers);
    
    void compute_optimality_measures(Problem& problem, Iterate& iterate);
    void compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution);
    
    /* use pointer to allow polymorphism */
    std::shared_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */

protected:
    /* QP subproblems */
    // optimality QP
    SubproblemSolution compute_qp_step(Problem& problem, Iterate& current_iterate, double trust_region_radius = INFINITY);
    double compute_qp_predicted_reduction(Iterate& current_iterate, SubproblemSolution& solution, double step_length);
    virtual std::vector<Range> generate_variables_bounds(Problem& problem, Iterate& current_iterate, double trust_region_radius);
    // feasibility QP
    SubproblemSolution compute_feasibility_qp_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius);
    void compute_linear_feasibility_objective(Iterate& current_iterate, ConstraintPartition& constraint_partition);
    std::vector<double> generate_feasibility_multipliers(Problem& problem, std::vector<double>& current_constraint_multipliers, ConstraintPartition& constraint_partition);
    std::vector<Range> generate_feasibility_bounds(Problem& problem, std::vector<double>& current_constraints, ConstraintPartition& constraint_partition);
    //virtual void evaluate_optimality_iterate(Problem& problem, Iterate& current_iterate) = 0;
    //virtual void evaluate_feasibility_iterate(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution) = 0;
    
    /* LP subproblems */
    SubproblemSolution compute_lp_step(Problem& problem, Iterate& current_iterate, double trust_region_radius);
    double compute_lp_predicted_reduction(SubproblemSolution& solution, double step_length);
    SubproblemSolution compute_feasibility_lp_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius);
};

#endif // ACTIVESETMETHOD_H
