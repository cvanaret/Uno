// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "options/Options.hpp"

namespace uno {
   ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Options& options):
         primal_tolerance(options.get_double("primal_tolerance")),
         dual_tolerance(options.get_double("dual_tolerance")),
         loose_primal_tolerance(options.get_double("loose_primal_tolerance")),
         loose_dual_tolerance(options.get_double("loose_dual_tolerance")),
         loose_tolerance_iteration_threshold(options.get_unsigned_int("loose_tolerance_iteration_threshold")),
         diverging_iterate_threshold(options.get_double("diverging_iterate_threshold")),
         unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")) {
   }

   ConstraintRelaxationStrategy::~ConstraintRelaxationStrategy() = default;

   size_t ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
      return this->number_subproblems_solved;
   }
} // namespace