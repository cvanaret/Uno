#ifndef SLP_H
#define SLP_H

#include "Subproblem.hpp"
#include "QPSolver.hpp"

class SLP : public Subproblem {
public:
   SLP(const Problem& problem, std::string QP_solver_name, bool use_trust_region, bool scale_residuals);

   void generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;

   Direction compute_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   int get_hessian_evaluation_count() const override;

private:
   /* use references to allow polymorphism */
   std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */

   static double compute_predicted_reduction_(Direction& direction, double step_length);
};

#endif // SLP_H
