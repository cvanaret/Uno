#ifndef SQP_H
#define SQP_H

#include "ActiveSetMethod.hpp"
#include "QPSolver.hpp"
#include "HessianEvaluation.hpp"

/*! \class QPApproximation
 * \brief QP local approximation
 *
 *  Quadratic approximation
 */
class SQP : public ActiveSetMethod {
public:
    /*!
     *  Constructor
     * 
     * \param solver: solver that solves the subproblem
     */
    SQP(Problem& problem, std::string QP_solver_name, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals);
    bool phase_1_required(SubproblemSolution& solution) override;

    /* use references to allow polymorphism */
    std::shared_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */

protected:
    void evaluate_optimality_iterate(Problem& problem, Iterate& current_iterate);
    // phase 1
    virtual void evaluate_feasibility_iterate(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution);
    std::vector<double> generate_feasibility_multipliers(Problem& problem, std::vector<double>& current_constraint_multipliers, ConstraintPartition& constraint_partition);
    std::vector<Range> generate_feasibility_bounds(Problem& problem, std::vector<double>& current_constraints, ConstraintPartition& constraint_partition);
    void set_feasibility_objective_(Iterate& current_iterate, ConstraintPartition& constraint_partition);
    // call subproblem solver
    virtual SubproblemSolution solve_optimality_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, Iterate& current_iterate, std::vector<double>& d0);
    virtual SubproblemSolution solve_feasibility_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, Iterate& current_iterate, std::vector<double>& d0);
    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length);
};

#endif // SQP_H
