// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BACKTRACKINGLINESEARCH_H
#define UNO_BACKTRACKINGLINESEARCH_H

#include "GlobalizationMechanism.hpp"

class BacktrackingLineSearch : public GlobalizationMechanism {
public:
   BacktrackingLineSearch(Statistics& statistics, ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);

   void initialize(Iterate& initial_iterate) override;
   [[nodiscard]] Iterate compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate) override;

private:
   const double backtracking_ratio;
   const double minimum_step_length;
   const bool scale_duals_with_step_length;
   size_t total_number_iterations{0}; /*!< Total number of iterations (optimality and feasibility) */

   [[nodiscard]] Iterate backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate, const Direction& direction,
      WarmstartInformation& warmstart_information);
   [[nodiscard]] Iterate assemble_trial_iterate(const Model& model, Iterate& current_iterate, const Direction& direction,
         double primal_dual_step_length) const;
   [[nodiscard]] double decrease_step_length(double step_length) const;
   static void check_unboundedness(const Direction& direction);
   void set_statistics(Statistics& statistics, const Direction& direction, double primal_dual_step_length) const;
   static void print_iteration(size_t number_iterations, double primal_dual_step_length);
};

#endif // UNO_BACKTRACKINGLINESEARCH_H
