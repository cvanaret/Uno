#ifndef L1RELAXATION_H
#define L1RELAXATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/strategy/GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

struct l1RelaxationParameters {
   bool fixed_parameter;
   double decrease_factor;
   double epsilon1;
   double epsilon2;
   double small_threshold;
};

class l1Relaxation : public ConstraintRelaxationStrategy {
public:
   l1Relaxation(Problem& problem, const Options& options);
   void initialize(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& first_iterate) override;

   // direction computation
   void create_current_subproblem(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, double trust_region_radius) override;
   Direction compute_feasible_direction(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate) override;
   Direction compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) override;
   Direction solve_feasibility_problem(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
         const Direction& phase_2_direction) override;

   bool is_acceptable(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, PredictedReductionModel& predicted_reduction_model, double step_length) override;

protected:
   const std::unique_ptr<GlobalizationStrategy> globalization_strategy;
   double penalty_parameter;
   const l1RelaxationParameters parameters;
   // preallocated temporary multipliers
   std::vector<double> constraint_multipliers;

   static double compute_predicted_reduction(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, const Direction& direction,
         PredictedReductionModel& predicted_reduction_model, double step_length);
   Direction solve_subproblem(Statistics& statistics, const Problem& problem, Iterate& current_iterate, double current_penalty_parameter);
   Direction resolve_subproblem(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
         double current_penalty_parameter);
   Direction solve_with_steering_rule(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate);
   double compute_linearized_constraint_residual(std::vector<double>& direction) const;
   double compute_error(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, const Multipliers& multipliers_displacements,
         double current_penalty_parameter);
   static void set_multipliers(const Problem& problem, const Iterate& current_iterate, std::vector<double>& constraints_multipliers);
};

#endif //L1RELAXATION_H