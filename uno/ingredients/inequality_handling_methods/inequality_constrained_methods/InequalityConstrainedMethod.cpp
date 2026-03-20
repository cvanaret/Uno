// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "InequalityConstrainedMethod.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/constraint_relaxation_strategies/relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/EvaluationCache.hpp"

namespace uno {
   void InequalityConstrainedMethod::check_problem(const OptimizationProblem& /*problem*/, bool /*uses_trust_region*/) {
   }

   void InequalityConstrainedMethod::initialize_statistics(Statistics& /*statistics*/) {
      // do nothing
   }

   std::unique_ptr<OptimizationProblem> InequalityConstrainedMethod::reformulate(const OptimizationProblem& problem,
         Parameterization& /*parameterization*/) {
      return problem.clone(); // the problem is not reformulated
   }

   bool InequalityConstrainedMethod::update_parameterization(Statistics& /*statistics*/, const OptimizationProblem& /*problem*/,
         const Iterate& /*current_iterate*/, Parameterization& /*parameterization*/) {
      return false;
   }

   void InequalityConstrainedMethod::initialize_feasibility_problem(Iterate& /*current_iterate*/) {
      // do nothing
   }

   void InequalityConstrainedMethod::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate,
         Evaluations& /*evaluations*/) {
      // c(x) - p + n = 0
      // TODO set (one of) the elastic variables to +/- the value of the constraint if violated
      problem.set_elastic_variable_values(current_iterate, [&](Iterate& iterate, size_t /*j*/, size_t elastic_index, double /*jacobian_coefficient*/) {
         iterate.primals[elastic_index] = 0.;
         iterate.multipliers.lower_bounds[elastic_index] = 1.;
         iterate.multipliers.upper_bounds[elastic_index] = 0.;
      });
   }

   double InequalityConstrainedMethod::proximal_coefficient() const {
      return 0.;
   }

   // compute dual *displacements*
   // because of the way we form LPs/QPs, we get the new *multipliers* back from the solver. To get the dual displacements/direction,
   // we need to subtract the current multipliers
   void InequalityConstrainedMethod::compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers) {
      view(direction_multipliers.constraints, 0, current_multipliers.constraints.size()) -= current_multipliers.constraints;
      view(direction_multipliers.lower_bounds, 0, current_multipliers.lower_bounds.size()) -= current_multipliers.lower_bounds;
      view(direction_multipliers.upper_bounds, 0, current_multipliers.upper_bounds.size()) -= current_multipliers.upper_bounds;
   }

   std::string InequalityConstrainedMethod::get_name() const {
      return "inequality-constrained method";
   }
} // namespace