#ifndef SLP_H
#define SLP_H

#include "ActiveSetMethod.hpp"

class SLP : public ActiveSetMethod {
public:
   SLP(Problem& problem, std::string QP_solver_name, bool use_trust_region, bool scale_residuals);

   void generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multipliers(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;

   std::vector<Direction> compute_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   std::vector<Direction>
   restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) override;

   /* use references to allow polymorphism */
   std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */

private:
   void evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate);
   void evaluate_feasibility_iterate_(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction);
};

#endif // SLP_H
