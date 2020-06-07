#ifndef SLPEQP_H
#define SLPEQP_H

#include <memory>
#include "ActiveSetMethod.hpp"
#include "QPSolver.hpp"
#include "LinearSolver.hpp"

class SLPEQP: public ActiveSetMethod {
public:
    SLPEQP(Problem& problem, std::string LP_solver_name, std::string linear_solver_name, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals);
    
    Direction compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius=INFINITY) override;
    Direction restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius=INFINITY) override;
    
    /* use pointers to allow polymorphism */
    std::shared_ptr<QPSolver> lp_solver;
    std::shared_ptr<LinearSolver> linear_solver; /*!< Solver that solves the subproblem */
    std::shared_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */

private:
    Direction solve_eqp_(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius);
    //void fix_active_constraints_(Problem& problem, ActiveSet& active_set, std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds);
    //virtual SubproblemSolution compute_optimality_eqp_step_() = 0;
    //virtual SubproblemSolution compute_feasibility_eqp_step_() = 0;
};

//class SLPEQP_TR : public SLPEQP {
//public:
//    SLPEQP_TR(Problem& problem, std::string LP_solver_name, std::string QP_solver_name, std::string hessian_evaluation_method, bool scale_residuals);
//    
//    /* use a pointer to allow polymorphism */
//    std::shared_ptr<QPSolver> qp_solver; /*!< Solver that solves the subproblem */
//};
//
//class SLPEQP_l2 : public SLPEQP {
//public:
//    SLPEQP_l2(Problem& problem, std::string LP_solver_name, std::string linear_solver_name, std::string hessian_evaluation_method, bool scale_residuals);
//};

#endif // SLPEQP_H
