#ifndef Sl1QP_H
#define Sl1QP_H

#include "Subproblem.hpp"
#include "QPSolver.hpp"
#include "HessianEvaluation.hpp"

struct Sl1QPParameters {
    double tau;
    double epsilon1;
    double epsilon2;
};

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
    Sl1QP(Problem& problem, std::string QP_solver, std::string hessian_evaluation_method);

    Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool use_trust_region) override;

    SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, double trust_region_radius = INFINITY) override;
    SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius = INFINITY) override;
    
    void compute_optimality_measures(Problem& problem, Iterate& iterate) override;
    void compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution);
    
    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) override;
    bool phase_1_required(SubproblemSolution& solution) override;
    double compute_complementarity_error(Problem& problem, Iterate& iterate, Multipliers& multipliers) override;

    /* use pointers to allow polymorphism */
    std::shared_ptr<QPSolver> solver; /*!< Subproblem solver */
    std::shared_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */
    double penalty_parameter;
    Sl1QPParameters parameters;

private:
    /* problem reformulation */
    // constraints l <= c(x) = u are reformulated as c(x) - p + n
    std::map<int, int> positive_part_variables; // p
    std::map<int, int> negative_part_variables; // n

    std::vector<Range> generate_variables_bounds(Problem& problem, Iterate& current_iterate, double trust_region_radius);
    SubproblemSolution solve_subproblem(Problem& problem, Iterate& current_iterate, double trust_region_radius, double penalty_parameter);
    double compute_linearized_constraint_residual(Problem& problem, std::vector<double>& x);
    double compute_error(Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter);
    void recover_active_set(Problem& problem, SubproblemSolution& solution, std::vector<Range>& variables_bounds);
};

#endif // Sl1QP_H
