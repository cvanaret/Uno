#ifndef SQP_H
#define SQP_H

#include "Subproblem.hpp"
#include "HessianEvaluation.hpp"
#include "QPSolver.hpp"

class SQP : public Subproblem {
public:
   SQP(const Problem& problem, size_t number_variables, size_t number_constraints, const std::string& QP_solver_name,
         const std::string& hessian_evaluation_method, bool use_trust_region);

   void generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;
   void set_initial_point(const std::vector<double>& point) override;

   Direction compute_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   double compute_predicted_reduction(const Direction& direction, double step_length) const override;
   int get_hessian_evaluation_count() const override;

protected:
   /* use references to allow polymorphism */
   std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */
   std::unique_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */
   std::vector<double> initial_point;
};

#endif // SQP_H
