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
   [[nodiscard]] Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, const std::vector<double>& initial_point,
         WarmstartInformation& warmstart_information) override;
   void switch_to_feasibility_problem(Iterate& current_iterate, WarmstartInformation& warmstart_information) override;

   // trial iterate acceptance
   void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length) override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) override;

   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] size_t get_number_subproblems_solved() const override;

protected:
   const l1RelaxedProblem feasibility_problem;
   l1RelaxedProblem l1_relaxed_problem;
   std::unique_ptr<Subproblem> subproblem;
   const std::unique_ptr<GlobalizationStrategy> globalization_strategy;
   double penalty_parameter;
   const double tolerance;
   const l1RelaxationParameters parameters;
   const double small_duals_threshold;
   // preallocated temporary multipliers
   Multipliers trial_multipliers;

   Direction solve_sequence_of_relaxed_subproblems(Statistics& statistics, Iterate& current_iterate, WarmstartInformation& warmstart_information);
   Direction solve_l1_relaxed_problem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter,
         const WarmstartInformation& warmstart_information);
   Direction solve_subproblem(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
         const WarmstartInformation& warmstart_information);

   // functions that decrease the penalty parameter to enforce particular conditions
   void decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction);
   double compute_infeasible_dual_error(Iterate& current_iterate);
   [[nodiscard]] Direction enforce_linearized_residual_sufficient_decrease(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double linearized_residual, double residual_lowest_violation, WarmstartInformation& warmstart_information);
   [[nodiscard]] bool linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual,
         double residual_lowest_violation) const;
   [[nodiscard]] Direction enforce_descent_direction_for_l1_merit(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         const Direction& direction_lowest_violation, WarmstartInformation& warmstart_information);
   [[nodiscard]] bool is_descent_direction_for_l1_merit_function(const Iterate& current_iterate, const Direction& direction,
         const Direction& direction_lowest_violation) const;

   void set_progress_measures(Iterate& iterate) const;
   [[nodiscard]] ProgressMeasures compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction,
         double step_length);

   [[nodiscard]] double compute_complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers) const override;

   void add_statistics(Statistics& statistics, const Iterate& trial_iterate) const;
   void check_exact_relaxation(Iterate& iterate) const;
};

#endif //UNO_L1RELAXATION_H