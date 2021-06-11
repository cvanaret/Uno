#ifndef SQP_H
#define SQP_H

#include "ActiveSetMethod.hpp"
#include "HessianEvaluation.hpp"

class SQP : public ActiveSetMethod {
public:
   SQP(const Problem& problem, const std::string& QP_solver_name, const std::string& hessian_evaluation_method, bool use_trust_region, bool scale_residuals);

   void generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multipliers(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;

   Direction compute_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   Direction restore_feasibility(const Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) override;
   int get_hessian_evaluation_count() override;

protected:
   /* use references to allow polymorphism */
   std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */
   std::unique_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */

   Direction compute_l1qp_step_(const Problem& problem, Iterate& current_iterate, ConstraintPartition& constraint_partition,
         std::vector<double>& initial_point, double trust_region_radius);
   void evaluate_feasibility_iterate_(const Problem& problem, Iterate& current_iterate, ConstraintPartition& constraint_partition);
   double compute_predicted_reduction_(Iterate& current_iterate, Direction& direction, double step_length) const;
};

#endif // SQP_H
