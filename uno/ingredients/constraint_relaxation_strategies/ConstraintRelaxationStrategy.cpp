// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/Options.hpp"

namespace uno {
   ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Options& options):
         residual_norm(norm_from_string(options.get_string("residual_norm"))),
         primal_tolerance(options.get_double("primal_tolerance")),
         dual_tolerance(options.get_double("dual_tolerance")),
         loose_primal_tolerance(options.get_double("loose_primal_tolerance")),
         loose_dual_tolerance(options.get_double("loose_dual_tolerance")),
         loose_tolerance_iteration_threshold(options.get_unsigned_int("loose_tolerance_iteration_threshold")),
         diverging_iterate_threshold(options.get_double("diverging_iterate_threshold")),
         unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")),
         acceptable_tolerance(1e-6),                    // scaled overall error [1e-6]
         acceptable_dual_infeasibility_tolerance(1e10), // unscaled [1e10]
         acceptable_primal_infeasibility_tolerance(1e-2), // unscaled [1e-2]
         acceptable_complementarity_tolerance(1e-2),    // unscaled [1e-2]
         acceptable_objective_change_tolerance(1e20),   // relative [1e20]
         acceptable_iteration_threshold(15) {          // [15], 0 disables the heuristic{
   }

   ConstraintRelaxationStrategy::~ConstraintRelaxationStrategy() = default;

   size_t ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
      return this->number_subproblems_solved;
   }
} // namespace