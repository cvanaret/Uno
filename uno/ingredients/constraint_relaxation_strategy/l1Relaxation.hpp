// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXATION_H
#define UNO_L1RELAXATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategy.hpp"
#include "reformulation/OptimalityProblem.hpp"
#include "reformulation/l1RelaxedProblem.hpp"
#include "tools/Options.hpp"

struct l1RelaxationParameters {
   bool fixed_parameter;
   double decrease_factor;
   double epsilon1;
   double epsilon2;
   double residual_small_threshold;
};

class l1Relaxation : public ConstraintRelaxationStrategy {
public:
   l1Relaxation(Statistics& statistics, const Model& model, const Options& options);
   void initialize(Iterate& initial_iterate) override;

   void set_trust_region_radius(double trust_region_radius) override;

   // direction computation
   [[nodiscard]] Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
         WarmstartInformation& warmstart_information) override;
   [[nodiscard]] Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         WarmstartInformation& warmstart_information) override;
   [[nodiscard]] Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         const std::vector<double>& initial_point, WarmstartInformation& warmstart_information) override;

   // trial iterate acceptance
   void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length) override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) override;

   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] size_t get_number_subproblems_solved() const override;

protected:
   const l1RelaxedProblem feasibility_problem;
   l1RelaxedProblem relaxed_problem;
   std::unique_ptr<Subproblem> subproblem;
   const std::unique_ptr<GlobalizationStrategy> globalization_strategy;
   double penalty_parameter;
   const l1RelaxationParameters parameters;
   const double small_duals_threshold{1e-8};
   const double l1_constraint_violation_coefficient;
   // preallocated temporary multipliers
   Multipliers trial_multipliers;

   Direction solve_subproblem(Statistics& statistics, Iterate& current_iterate, const NonlinearProblem& problem,
         const WarmstartInformation& warmstart_information);
   Direction solve_relaxed_problem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter,
         const WarmstartInformation& warmstart_information);
   Direction solve_with_steering_rule(Statistics& statistics, Iterate& current_iterate, WarmstartInformation& warmstart_information);
   void decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction);
   [[nodiscard]] bool linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual, double residual_lowest_violation) const;
   [[nodiscard]] bool is_descent_direction_for_l1_merit_function(const Iterate& current_iterate, const Direction& direction, const Direction& direction_lowest_violation) const;
   double compute_dual_error(Iterate& current_iterate);
   void check_exact_relaxation(Iterate& iterate) const;

   // progress measures and their local models
   void set_infeasibility_measure(Iterate& iterate);
   [[nodiscard]] double generate_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Direction& direction,
         double step_length) const;
   void set_optimality_measure(Iterate& iterate);
   [[nodiscard]] std::function<double (double)> generate_predicted_optimality_reduction_model(const Iterate& current_iterate,
         const Direction& direction, double step_length) const;
};

#endif //UNO_L1RELAXATION_H