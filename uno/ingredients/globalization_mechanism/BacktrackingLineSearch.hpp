// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BACKTRACKINGLINESEARCH_H
#define UNO_BACKTRACKINGLINESEARCH_H

#include "GlobalizationMechanism.hpp"

namespace uno {
   class BacktrackingLineSearch : public GlobalizationMechanism {
   public:
      BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);

      void initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) override;
      void compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate) override;

   private:
      const double backtracking_ratio;
      const double minimum_step_length;
      const bool scale_duals_with_step_length;
      size_t total_number_iterations{0}; /*!< Total number of iterations (optimality and feasibility) */

      void backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
            WarmstartInformation& warmstart_information);
      [[nodiscard]] double decrease_step_length(double step_length) const;
      static void check_unboundedness(const Direction& direction);
      void set_statistics(Statistics& statistics) const;
      void set_statistics(Statistics& statistics, const Iterate& trial_iterate, const Direction& direction, double primal_dual_step_length) const;
      static void print_iteration(size_t number_iterations, double primal_dual_step_length);
   };
} // namespace

#endif // UNO_BACKTRACKINGLINESEARCH_H
