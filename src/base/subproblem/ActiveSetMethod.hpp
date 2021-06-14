#ifndef ACTIVESETMETHOD_H
#define ACTIVESETMETHOD_H

#include <memory>
#include "Subproblem.hpp"
#include "QPSolver.hpp"
#include "HessianEvaluation.hpp"

struct ElasticVariables {
   std::map<int, int> positive;
   std::map<int, int> negative;
   [[nodiscard]] size_t size() const { return this->positive.size() + this->negative.size(); }
};

class ActiveSetMethod : public Subproblem {
public:
   ActiveSetMethod(const Problem& problem, bool scale_residuals);

   Iterate evaluate_initial_point(const Problem& problem, const std::vector<double>& x, const Multipliers& multipliers) override;
   void compute_optimality_measures(const Problem& problem, Iterate& iterate) override;
   static void generate_elastic_variables_(const Problem& problem, ElasticVariables& elastic_variables);

protected:
   static void recover_l1qp_active_set_(const Problem& problem, Direction& direction, const ElasticVariables& elastic_variables);

   /* LP subproblems */
   Direction compute_lp_step_(const Problem& problem, QPSolver& solver, Iterate& current_iterate, double trust_region_radius);
   static double compute_lp_predicted_reduction_(Direction& direction, double step_length);
   Direction compute_l1lp_step_(const Problem& problem, QPSolver& solver, Iterate& current_iterate, Direction& phase_2_direction,
         double trust_region_radius);
};

#endif // ACTIVESETMETHOD_H
