// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "ingredients/subproblem/PredictedReductionModel.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Statistics.hpp"
#include "tools/Options.hpp"

class ConstraintRelaxationStrategy {
public:
   ConstraintRelaxationStrategy(bool penalty_parameter_control, const Options& options);
   virtual ~ConstraintRelaxationStrategy() = default;

   virtual void initialize(Statistics& statistics, Iterate& first_iterate) = 0;
   virtual void set_variable_bounds(const Iterate& current_iterate, double trust_region_radius) = 0;

   // direction computation
   virtual Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) = 0;
   virtual Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         const std::optional<std::vector<double>>& optional_phase_2_solution) = 0;
   virtual Direction compute_second_order_correction(Iterate& trial_iterate) = 0;

   // trial iterate acceptance
   [[nodiscard]] virtual PredictedReductionModel generate_predicted_reduction_model(const Direction& direction) const = 0;
   virtual bool is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         PredictedReductionModel& predicted_reduction_model, double step_length = 1.) = 0;
   virtual void register_accepted_iterate(Iterate& iterate) = 0;

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   [[nodiscard]] virtual size_t get_number_subproblems_solved() const = 0;

   const bool penalty_parameter_control;

protected:
   const Norm residual_norm;
   const double small_step_threshold;

   [[nodiscard]] virtual double compute_infeasibility_measure(Iterate& iterate) = 0;
   [[nodiscard]] bool is_small_step(const Direction& direction) const;
   void compute_nonlinear_residuals(const ReformulatedProblem& problem, Iterate& iterate) const;
   void recover_active_set(const Model& model, Direction& direction);
};

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H
