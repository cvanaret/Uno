#ifndef SLPEQP_H
#define SLPEQP_H

#include <memory>
#include "ActiveSetMethod.hpp"
#include "QPSolver.hpp"
#include "LinearSolver.hpp"

class SLPEQP : public ActiveSetMethod {
public:
   SLPEQP(Problem& problem, std::string LP_solver_name, std::string linear_solver_name, std::string hessian_evaluation_method,
         bool use_trust_region, bool scale_residuals);

   void generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multipliers(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;

   std::vector<Direction> compute_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   std::vector<Direction>
   restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) override;

private:
   /* use pointers to allow polymorphism */
   std::unique_ptr<QPSolver> lp_solver;
   std::unique_ptr<LinearSolver> linear_solver; /*!< Solver that solves the subproblem */
   std::unique_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */

   Direction solve_eqp_(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius);
   double compute_qp_predicted_reduction_(Iterate& current_iterate, Direction& direction, double step_length);
   //void fix_active_constraints_(Problem& problem, ActiveSet& active_set, std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds);
   //virtual SubproblemSolution compute_optimality_eqp_step_() = 0;
   //virtual SubproblemSolution compute_feasibility_eqp_step_() = 0;
};

//class SLPEQP_TR : public SLPEQP {
//public:
//    SLPEQP_TR(Problem& problem, std::string LP_solver_name, std::string QP_solver_name, std::string hessian_evaluation_method, bool scale_residuals);
//    
//    /* use a pointer to allow polymorphism */
//    std::unique_ptr<QPSolver> qp_solver; /*!< Solver that solves the subproblem */
//};
//
//class SLPEQP_l2 : public SLPEQP {
//public:
//    SLPEQP_l2(Problem& problem, std::string LP_solver_name, std::string linear_solver_name, std::string hessian_evaluation_method, bool scale_residuals);
//};

#endif // SLPEQP_H
