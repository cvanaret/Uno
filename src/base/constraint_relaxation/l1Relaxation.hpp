#ifndef FEASIBILITYRESTORATION_H
#define FEASIBILITYRESTORATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "ActiveSetMethod.hpp"

struct RelaxationParameters {
   double tau;
   double epsilon1;
   double epsilon2;
};

class l1Relaxation: public ConstraintRelaxationStrategy {
public:
   explicit l1Relaxation(Problem& problem, Subproblem& subproblem);
   std::vector<Direction> compute_feasible_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   std::optional<Iterate> check_acceptance(Statistics& statistics, Problem& problem, Iterate& current_iterate, Direction& direction, double
   step_length) override;
   double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) override;

protected:
   double penalty_parameter;
   RelaxationParameters parameters;
   /* problem reformulation with elastic variables. Constraints l <= c(x) = u are reformulated as c(x) - p + n */
   ElasticVariables elastic_variables;

   void preprocess_subproblem();
   std::vector<Direction> compute_byrd_steering_rule(Problem& problem, Iterate& current_iterate, double trust_region_radius);
   double compute_linearized_constraint_residual(std::vector<double>& direction);
   double compute_error(Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter);
   void postprocess_direction(const Problem& problem, Direction& direction);
   double compute_complementarity_error(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) const;
};

#endif //FEASIBILITYRESTORATION_H