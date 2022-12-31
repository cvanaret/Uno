// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include "ingredients/globalization_strategy/PredictedReductionModel.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Statistics.hpp"
#include "tools/Options.hpp"

class ConstraintRelaxationStrategy {
public:
   ConstraintRelaxationStrategy(const Model& model, const Options& options);
   virtual ~ConstraintRelaxationStrategy() = default;

   virtual void initialize(Statistics& statistics, Iterate& first_iterate) = 0;
   virtual void set_trust_region_radius(double trust_region_radius) = 0;

   // direction computation
   [[nodiscard]] virtual Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) = 0;
   [[nodiscard]] virtual Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate) = 0;
   [[nodiscard]] virtual Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         const std::vector<double>& initial_point) = 0;
   [[nodiscard]] virtual Direction compute_second_order_correction(Iterate& trial_iterate) = 0;

   // trial iterate acceptance
   virtual void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction) = 0;
   [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) = 0;
   virtual void register_accepted_iterate(Iterate& iterate) = 0;

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   [[nodiscard]] virtual size_t get_number_subproblems_solved() const = 0;

protected:
   const Model& original_model;
   const Norm residual_norm;
   const double small_step_threshold;

   [[nodiscard]] static double compute_linearized_constraint_violation(const Model& model, const Iterate& current_iterate, const Direction& direction,
         double step_length);
   [[nodiscard]] bool is_small_step(const Direction& direction) const;
   static void evaluate_lagrangian_gradient(Iterate& iterate, const std::vector<double>& constraint_multipliers,
         const std::vector<double>& lower_bound_multipliers, const std::vector<double>& upper_bound_multipliers);
   static void evaluate_reformulation_functions(const NonlinearProblem& problem, Iterate& iterate);
   void compute_primal_dual_errors(const NonlinearProblem& problem, Iterate& iterate) const;
};

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H
