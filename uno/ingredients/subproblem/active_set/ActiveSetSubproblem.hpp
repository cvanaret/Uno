// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ACTIVESETSUBPROBLEM_H
#define UNO_ACTIVESETSUBPROBLEM_H

#include "ingredients/subproblem/Subproblem.hpp"

class ActiveSetSubproblem : public Subproblem {
public:
   ActiveSetSubproblem(size_t max_number_variables, size_t max_number_constraints);
   ~ActiveSetSubproblem() override = default;

   void generate_initial_iterate(const NonlinearProblem& problem, Iterate& initial_iterate) override;
   void set_initial_point(const std::vector<double>& initial_point) override;
   void initialize_feasibility_problem() override;
   void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
   void exit_feasibility_problem(const NonlinearProblem& problem, Iterate& trial_iterate) override;

   void set_auxiliary_measure(const NonlinearProblem& problem, Iterate& iterate) override;
   [[nodiscard]] double generate_predicted_auxiliary_reduction_model(const NonlinearProblem&, const Iterate&, const Direction&, double) const override;

   void postprocess_iterate(const NonlinearProblem& model, Iterate& iterate) override;

protected:
   std::vector<double> initial_point{};
   std::vector<Interval> direction_bounds{};
   std::vector<Interval> linearized_constraint_bounds{};

   void set_direction_bounds(const NonlinearProblem& problem, const Iterate& current_iterate);
   void set_linearized_constraint_bounds(const NonlinearProblem& problem, const std::vector<double>& current_constraints);
   static void compute_dual_displacements(const NonlinearProblem& problem, const Iterate& current_iterate, Direction& direction);
};

#endif // UNO_ACTIVESETSUBPROBLEM_H