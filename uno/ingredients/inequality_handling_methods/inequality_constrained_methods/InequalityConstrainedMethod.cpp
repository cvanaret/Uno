// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityConstrainedMethod.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/constraint_relaxation_strategies/relaxed_problems/l1RelaxedProblem.hpp"

namespace uno {
   void InequalityConstrainedMethod::check_problem(const OptimizationProblem& /*problem*/, bool /*uses_trust_region*/) {
   }

   void InequalityConstrainedMethod::initialize_statistics(Statistics& /*statistics*/) {
      // do nothing
   }

   std::unique_ptr<OptimizationProblem> InequalityConstrainedMethod::reformulate(const OptimizationProblem& problem,
         Parameterization& /*parameterization*/) {
      return problem.clone(); // the problem is not reformulated, simply make a copy
   }

   bool InequalityConstrainedMethod::update_parameterization(Statistics& /*statistics*/, const OptimizationProblem& /*problem*/,
         const Iterate& /*current_iterate*/, Parameterization& /*parameterization*/) {
      // the parameterization is not updated
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

   std::string InequalityConstrainedMethod::get_name() const {
      return "inequality-constrained method";
   }
} // namespace