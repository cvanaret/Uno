#ifndef SQP_H
#define SQP_H

#include "ActiveSetMethod.hpp"
#include "HessianEvaluation.hpp"

class SQP : public ActiveSetMethod {
public:
   SQP(Problem& problem, std::string QP_solver_name, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals);

   void evaluate_current_iterate(const Problem& problem, const Iterate& current_iterate) override;

   std::vector<Direction> compute_directions(Problem& problem, Iterate& current_iterate, double objective_multiplier,
         double trust_region_radius) override;
   std::vector<Direction> restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction,
         double trust_region_radius) override;
   double compute_qp_predicted_reduction_(Problem& problem, Iterate& current_iterate, Direction& direction, double step_length);

   /* use references to allow polymorphism */
   std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */
   std::unique_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */

protected:
   Direction compute_qp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, double trust_region_radius);
   void evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate) const;
   Direction compute_l1qp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, ConstraintPartition& constraint_partition,
         std::vector<double>& initial_solution, double trust_region_radius);
   void evaluate_feasibility_iterate_(Problem& problem, Iterate& current_iterate, ConstraintPartition& constraint_partition);
};

#endif // SQP_H
