#ifndef L1RELAXATION_H
#define L1RELAXATION_H

#include <GlobalizationStrategy.hpp>
#include "ConstraintRelaxationStrategy.hpp"

struct RelaxationParameters {
   double tau;
   double epsilon1;
   double epsilon2;
};

struct ElasticVariables {
   std::map<size_t, size_t> positive;
   std::map<size_t, size_t> negative;
   [[nodiscard]] size_t size() const { return this->positive.size() + this->negative.size(); }
};

class l1Relaxation: public ConstraintRelaxationStrategy {
public:
   explicit l1Relaxation(Problem& problem, Subproblem& subproblem, const std::map<std::string, std::string>& options);
   Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;

   // direction computation
   void generate_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   Direction compute_feasible_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   Direction solve_feasibility_problem(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Direction& direction) override;

   bool is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) override;
   double compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, const Direction& direction, double step_length) override;

protected:
   std::unique_ptr<GlobalizationStrategy> globalization_strategy;
   double penalty_parameter;
   RelaxationParameters parameters;
   /* problem reformulation with elastic variables. Constraints l <= c(x) = u are reformulated as l <= c(x) - p + n <= u */
   ElasticVariables elastic_variables;

   Direction compute_byrd_steering_rule(Statistics& statistics, const Problem& problem, Iterate& current_iterate);
   void generate_elastic_variables_(const Problem& problem);
   double compute_linearized_constraint_residual(std::vector<double>& direction);
   double compute_error(const Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter);
   void postprocess_direction(const Problem& problem, Direction& direction);
   [[nodiscard]] double compute_complementarity_error(const Problem& problem, const Iterate& iterate, const Multipliers& multipliers) const;
   void recover_l1qp_active_set_(const Problem& problem, const Direction& direction);
};

#endif //L1RELAXATION_H