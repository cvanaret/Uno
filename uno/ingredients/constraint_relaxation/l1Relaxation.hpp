#ifndef L1RELAXATION_H
#define L1RELAXATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/strategy/GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

struct RelaxationParameters {
   double tau;
   double epsilon1;
   double epsilon2;
};

class l1Relaxation : public ConstraintRelaxationStrategy {
public:
   l1Relaxation(Problem& problem, const Options& options, bool use_trust_region);
   Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;

   // direction computation
   void generate_subproblem(const Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   Direction compute_feasible_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   Direction solve_feasibility_problem(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Direction& direction) override;

   bool is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) override;
   double compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, const Direction& direction, double step_length) override;

protected:
   const std::unique_ptr <GlobalizationStrategy> globalization_strategy;
   double penalty_parameter;
   const RelaxationParameters parameters;
   /* problem reformulation with elastic variables. Constraints l <= c(x) = u are reformulated as l <= c(x) - p + n <= u */
   ElasticVariables elastic_variables;

   Direction solve_subproblem(Statistics& statistics, const Problem& problem, Iterate& current_iterate);
   Direction solve_subproblem(Statistics& statistics, const Problem& problem, Iterate& current_iterate, double objective_multiplier);
   Direction compute_byrd_steering_rule(Statistics& statistics, const Problem& problem, Iterate& current_iterate);
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier);
   static size_t count_elastic_variables(const Problem& problem);
   double compute_linearized_constraint_residual(std::vector<double>& direction);
   double compute_error(const Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter) const;
   void remove_elastic_variables(const Problem& problem, Direction& direction);
   void recover_l1qp_active_set_(const Problem& problem, const Direction& direction);
};

#endif //L1RELAXATION_H