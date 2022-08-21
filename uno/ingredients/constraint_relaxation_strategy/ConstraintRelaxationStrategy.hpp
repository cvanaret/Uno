// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "ingredients/subproblem/PredictedOptimalityReductionModel.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Statistics.hpp"
#include "tools/Options.hpp"

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
   virtual bool is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         PredictedOptimalityReductionModel& predicted_optimality_reduction_model, double step_length) = 0;
   [[nodiscard]] virtual PredictedOptimalityReductionModel generate_predicted_optimality_reduction_model(const Direction& direction) const = 0;
   virtual void register_accepted_iterate(Iterate& iterate) = 0;

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   [[nodiscard]] virtual size_t get_number_subproblems_solved() const = 0;

protected:
   const Model& original_model;
   const Norm residual_norm;
   const double small_step_threshold;

   virtual void set_infeasibility_measure(Iterate& iterate) = 0;
   [[nodiscard]] static double compute_predicted_infeasibility_reduction(const Model& model, const Iterate& current_iterate, const Direction& direction,
         double step_length);
   [[nodiscard]] bool is_small_step(const Direction& direction) const;
   void compute_nonlinear_residuals(const NonlinearProblem& problem, Iterate& iterate) const;
};

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H
