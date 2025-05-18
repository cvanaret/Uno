// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityConstrainedMethod.hpp"
#include "optimization/Iterate.hpp"
#include "linear_algebra/Vector.hpp"
#include "ingredients/constraint_relaxation_strategies/l1RelaxedProblem.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   InequalityConstrainedMethod::InequalityConstrainedMethod(): InequalityHandlingMethod() {
   }

   std::pair<size_t, size_t> InequalityConstrainedMethod::get_dimensions(const OptimizationProblem& problem) const {
      return {problem.number_variables, problem.number_constraints};
   }

   void InequalityConstrainedMethod::initialize(const OptimizationProblem& problem, const HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy) {
      this->initial_point.resize(problem.number_variables);
      regularization_strategy.initialize_memory(problem, hessian_model);
   }

   void InequalityConstrainedMethod::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) {
   }

   void InequalityConstrainedMethod::set_initial_point(const Vector<double>& point) {
      // copy the point into the member
      this->initial_point = point;
   }

   void InequalityConstrainedMethod::initialize_feasibility_problem(const l1RelaxedProblem& /*problem*/, Iterate& /*current_iterate*/) {
      // do nothing
   }

   void InequalityConstrainedMethod::set_elastic_variable_values(const l1RelaxedProblem& problem, Vector<double>& current_primals,
         Multipliers& current_multipliers) {
      problem.set_elastic_variable_values([&](size_t /*j*/, size_t elastic_index, double /*jacobian_coefficient*/) {
         current_primals[elastic_index] = 0.;
         current_multipliers.lower_bounds[elastic_index] = 1.;
         current_multipliers.upper_bounds[elastic_index] = 0.;
      });
   }

   double InequalityConstrainedMethod::proximal_coefficient() const {
      return 0.;
   }

   void InequalityConstrainedMethod::exit_feasibility_problem(const OptimizationProblem& /*problem*/, Iterate& /*trial_iterate*/) {
      // do nothing
   }

   void InequalityConstrainedMethod::compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers) {
      // compute dual *displacements* (active-set methods usually compute the new duals, not the displacements)
      view(direction_multipliers.constraints, 0, current_multipliers.constraints.size()) -= current_multipliers.constraints;
      view(direction_multipliers.lower_bounds, 0, current_multipliers.lower_bounds.size()) -= current_multipliers.lower_bounds;
      view(direction_multipliers.upper_bounds, 0, current_multipliers.upper_bounds.size()) -= current_multipliers.upper_bounds;
   }

   // auxiliary measure is 0 in inequality-constrained methods
   double InequalityConstrainedMethod::compute_auxiliary_measure(const OptimizationProblem& /*first_reformulation*/, Iterate& /*iterate*/) {
      return 0.;
   }

   double InequalityConstrainedMethod::compute_predicted_auxiliary_reduction_model(const Model& /*model*/, const Iterate& /*current_iterate*/,
         const Vector<double>& /*primal_direction*/, double /*step_length*/) const {
      return 0.;
   }

   void InequalityConstrainedMethod::postprocess_iterate(const OptimizationProblem& /*problem*/, Vector<double>& /*primals*/, Multipliers& /*multipliers*/) {
   }
} // namespace
