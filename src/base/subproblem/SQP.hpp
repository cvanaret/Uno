#ifndef SQP_H
#define SQP_H

#include "Subproblem.hpp"
#include "HessianEvaluation.hpp"
#include "QPSolver.hpp"

class SQP : public Subproblem {
public:
   SQP(const Problem& problem, const std::string& QP_solver_name, const std::string& hessian_evaluation_method, bool use_trust_region, bool scale_residuals);

   void generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;

   Direction compute_direction(const Problem& problem, Iterate& current_iterate) override;
   int get_hessian_evaluation_count() override;

protected:
   /* use references to allow polymorphism */
   std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */
   std::unique_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */

   double compute_predicted_reduction_(Iterate& current_iterate, Direction& direction, double step_length) const;
};

#endif // SQP_H
