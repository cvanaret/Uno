#ifndef Sl1QP_H
#define Sl1QP_H

#include "Subproblem.hpp"
#include "QPSolver.hpp"
#include "HessianEvaluation.hpp"
#include "ActiveSetMethod.hpp"

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
class Sl1QP : public ActiveSetMethod {
public:
    /*!
     *  Constructor
     * 
     * \param solver: solver that solves the subproblem
     */
    Sl1QP(Problem& problem, std::string QP_solver, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals);
    
    SubproblemSolution compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius = INFINITY) override;
    SubproblemSolution restore_feasibility(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius = INFINITY) override;
    
    void evaluate_optimality_iterate(Problem& problem, Iterate& current_iterate, double penalty_parameter);
    
    double compute_complementarity_error(Problem& problem, Iterate& iterate, Multipliers& multipliers) override;

    /* use pointers to allow polymorphism */
    std::shared_ptr<QPSolver> solver; /*!< Subproblem solver */
    std::shared_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */
    double penalty_parameter;
    Sl1QPParameters parameters;
    int number_variables;

private:
    Sl1QP(Problem& problem, std::string QP_solver, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals, int number_variables);
    
    /* problem reformulation */
    // constraints l <= c(x) = u are reformulated as c(x) - p + n
    std::map<int, int> positive_part_variables; // p
    std::map<int, int> negative_part_variables; // n

    int count_additional_variables(Problem& problem);
    std::vector<Range> generate_variables_bounds(Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
    SubproblemSolution solve_l1qp_subproblem(Problem& problem, Iterate& current_iterate, double trust_region_radius, double penalty_parameter);
    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length);
    double compute_linearized_constraint_residual(std::vector<double>& x);
    double compute_error(Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter);
    void recover_active_set(Problem& problem, SubproblemSolution& solution, std::vector<Range>& variables_bounds);
};

#endif // Sl1QP_H
