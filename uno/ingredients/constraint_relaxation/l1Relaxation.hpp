#ifndef L1RELAXATION_H
#define L1RELAXATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/strategy/GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

struct l1RelaxationParameters {
   double initial_parameter;
   double decrease_factor;
   double epsilon1;
   double epsilon2;
};

class l1Relaxation : public ConstraintRelaxationStrategy {
public:
   l1Relaxation(Problem& problem, Subproblem& subproblem, const l1RelaxationParameters& parameters, const Options& options);
   void initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) override;

   // direction computation
   void create_current_subproblem(const Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   Direction compute_feasible_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   Direction solve_feasibility_problem(Statistics& statistics, const Problem& problem, Iterate& current_iterate, const Direction& phase_2_direction) override;

   bool is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         PredictedReductionModel& predicted_reduction_model, double step_length) override;
   double compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, const Direction& direction, PredictedReductionModel&
   predicted_reduction_model, double step_length) override;

   static size_t get_number_variables(const Problem& problem);

protected:
   const std::unique_ptr <GlobalizationStrategy> globalization_strategy;
   double penalty_parameter;
   const l1RelaxationParameters parameters;

   Direction solve_subproblem(Statistics& statistics, const Problem& problem, Iterate& current_iterate);
   Direction resolve_subproblem(Statistics& statistics, const Problem& problem, Iterate& current_iterate, double objective_multiplier);
   Direction solve_with_steering_rule(Statistics& statistics, const Problem& problem, Iterate& current_iterate);
   double compute_linearized_constraint_residual(std::vector<double>& direction) const;
   double compute_error(const Problem& problem, Iterate& iterate, Multipliers& multipliers, double current_penalty_parameter) const;
   void remove_elastic_variables(const Problem& problem, Direction& direction);
   void recover_active_set(const Problem& problem, const Direction& direction);
};

#endif //L1RELAXATION_H