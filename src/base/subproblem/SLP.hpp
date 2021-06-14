#ifndef SLP_H
#define SLP_H

#include "ActiveSetMethod.hpp"

class SLP : public ActiveSetMethod {
public:
   SLP(const Problem& problem, std::string QP_solver_name, bool use_trust_region, bool scale_residuals);

   void generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;

   Direction compute_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   //Direction restore_feasibility(const Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius)
   //override;
   int get_hessian_evaluation_count() override;

private:
   /* use references to allow polymorphism */
   std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */

   static double compute_predicted_reduction_(Iterate& current_iterate, Direction& direction, double step_length) ;
   void evaluate_optimality_iterate_(const Problem& problem, Iterate& current_iterate);
   void evaluate_feasibility_iterate_(const Problem& problem, Iterate& current_iterate, Direction& phase_2_direction);
};

#endif // SLP_H
