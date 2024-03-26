// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDMETHOD_H
#define UNO_INEQUALITYCONSTRAINEDMETHOD_H

#include "ingredients/subproblem/Subproblem.hpp"

class InequalityConstrainedMethod : public Subproblem {
public:
   InequalityConstrainedMethod(size_t max_number_variables, size_t max_number_constraints);
   ~InequalityConstrainedMethod() override = default;
   
   void initialize_statistics(Statistics& statistics, const Options& options) override;
   void set_initial_point(const std::vector<double>& point) override;
   void initialize_feasibility_problem() override;
   void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
   void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) override;

   void set_auxiliary_measure(const OptimizationProblem& problem, Iterate& iterate) override;
   [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const OptimizationProblem&, const Iterate&, const Direction&, double) const override;

   void postprocess_iterate(const OptimizationProblem& model, Iterate& iterate) override;

protected:
   std::vector<double> initial_point{};
   std::vector<Interval> direction_bounds{};
   std::vector<Interval> linearized_constraint_bounds{};

   void set_direction_bounds(const OptimizationProblem& problem, const Iterate& current_iterate);
   void set_linearized_constraint_bounds(const OptimizationProblem& problem, const std::vector<double>& current_constraints);
   static void compute_dual_displacements(const OptimizationProblem& problem, const Iterate& current_iterate, Direction& direction);
};

#endif // UNO_INEQUALITYCONSTRAINEDMETHOD_H
