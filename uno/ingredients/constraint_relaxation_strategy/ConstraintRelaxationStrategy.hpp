// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include <functional>
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Statistics.hpp"
#include "tools/Options.hpp"

struct PredictedReductionModel {
   std::function<double(double)> infeasibility;
   std::function<double(double)> scaled_optimality;
   std::function<double(double)> unscaled_optimality;
};

class ConstraintRelaxationStrategy {
public:
   ConstraintRelaxationStrategy(const Model& model, const Options& options);
   virtual ~ConstraintRelaxationStrategy() = default;

   virtual void initialize(Statistics& statistics, Iterate& first_iterate) = 0;
   virtual void set_variable_bounds(const Iterate& current_iterate, double trust_region_radius) = 0;

   // direction computation
   virtual Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) = 0;
   virtual Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate) = 0;
   virtual Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate, const std::vector<double>& initial_point) = 0;
   virtual Direction compute_second_order_correction(Iterate& trial_iterate) = 0;

   // trial iterate acceptance
   virtual void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction) = 0;
   virtual bool is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         PredictedReductionModel& predicted_reduction_model, double step_length) = 0;
   [[nodiscard]] virtual PredictedReductionModel generate_predicted_reduction_model(const Iterate& current_iterate, const Direction& direction) const = 0;
   virtual void register_accepted_iterate(Iterate& iterate) = 0;

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   [[nodiscard]] virtual size_t get_number_subproblems_solved() const = 0;

protected:
   const Model& original_model;
   const Norm residual_norm;
   const double small_step_threshold;

   static double compute_linearized_constraint_violation(const Model& model, const Iterate& current_iterate, const Direction& direction,
         double step_length);
   [[nodiscard]] bool is_small_step(const Direction& direction) const;
   static void evaluate_lagrangian_gradient(Iterate& iterate, const std::vector<double>& constraint_multipliers,
         const std::vector<double>& lower_bound_multipliers, const std::vector<double>& upper_bound_multipliers);
   static void evaluate_reformulation_functions(const NonlinearProblem& problem, Iterate& iterate);
   void compute_optimality_condition_residuals(const NonlinearProblem& problem, Iterate& iterate) const;
};

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H
