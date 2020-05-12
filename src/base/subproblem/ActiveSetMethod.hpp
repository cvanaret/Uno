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

    SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double trust_region_radius = INFINITY);
    SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius = INFINITY);
    void compute_optimality_measures(Problem& problem, Iterate& iterate);
    void compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution);
    virtual double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) = 0;
    std::vector<Range> generate_variables_bounds(std::vector<double>& current_x, std::vector<Range>& variables_bounds, double trust_region_radius);
    virtual bool phase_1_required(SubproblemSolution& solution) = 0;

    /* use pointer to allow polymorphism */
    std::shared_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */

protected:
    virtual void evaluate_optimality_iterate(Problem& problem, Iterate& current_iterate) = 0;
    // phase 1
    virtual void evaluate_feasibility_iterate(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution) = 0;
    std::vector<Range> generate_feasibility_bounds(Problem& problem, std::vector<double>& current_constraints, ConstraintPartition& constraint_partition);
    void compute_linear_feasibility_objective(Iterate& current_iterate, ConstraintPartition& constraint_partition);
    // call subproblem solver
    virtual SubproblemSolution solve_optimality_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, Iterate& current_iterate, std::vector<double>& d0) = 0;
    virtual SubproblemSolution solve_feasibility_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, Iterate& current_iterate, std::vector<double>& d0) = 0;
};

#endif // ACTIVESETMETHOD_H
